source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(survival)
  library(tibble)
})

# This is a sensitivity analysis, not a replacement candidate-definition rule.
# It removes the upstream hard LFC gate, models every QC-filtered gene available
# in the VST matrix, and applies one BH correction across the resulting Cox tests.
vst_mat <- read_required_rds(FILES$tcga_vst)
coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
deg <- read_csv(FILES$tcga_deg, show_col_types = FALSE)

required_deg_columns <- c("gene_id", "log2FoldChange", "log2FoldChange_apeglm", "lfcSE_apeglm")
missing_deg_columns <- setdiff(required_deg_columns, names(deg))
if (length(missing_deg_columns) > 0) {
  stop(
    "TCGA DESeq2 output is missing apeglm columns: ",
    paste(missing_deg_columns, collapse = ", "),
    ". Rerun analysis/04_deg_tcga.R."
  )
}

clinical_surv <- clinical |>
  transmute(
    patient_barcode = submitter_id,
    os_time = os_time,
    os_event = os_event,
    age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
    sex = factor(gender),
    stage = factor(normalize_stage(ajcc_pathologic_stage)),
    grade = factor(case_when(
      normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
      normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
      TRUE ~ NA_character_
    ))
  )

tumor_samples <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  filter(
    sample_barcode %in% colnames(vst_mat),
    !is.na(os_time), os_time > 0, !is.na(os_event),
    !is.na(age), !is.na(sex), !is.na(stage), !is.na(grade)
  ) |>
  mutate(
    sex = droplevels(sex),
    stage = droplevels(stage),
    grade = droplevels(grade)
  )

if (nrow(tumor_samples) < 100 || sum(tumor_samples$os_event == 1) < 25) {
  stop("Too few complete TCGA samples or events for the all-gene sensitivity analysis.")
}

fit_gene <- function(gene_id) {
  raw_expr <- vst_mat[gene_id, tumor_samples$sample_barcode]
  expr_sd <- stats::sd(raw_expr)
  if (!is.finite(expr_sd) || expr_sd == 0) {
    return(tibble(
      gene_id = gene_id,
      log_hr = NA_real_, hr = NA_real_, log_hr_ci_low = NA_real_,
      log_hr_ci_high = NA_real_, hr_ci_low = NA_real_, hr_ci_high = NA_real_,
      std_error = NA_real_, statistic = NA_real_, p_value = NA_real_,
      n = nrow(tumor_samples), events = sum(tumor_samples$os_event == 1),
      model_status = "zero_variance_expression", warning_text = NA_character_
    ))
  }

  dat <- tumor_samples |>
    transmute(
      os_time, os_event,
      expr = as.numeric(scale(raw_expr)),
      age, sex, stage, grade
    )

  warning_messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade, data = dat),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      gene_id = gene_id,
      log_hr = NA_real_, hr = NA_real_, log_hr_ci_low = NA_real_,
      log_hr_ci_high = NA_real_, hr_ci_low = NA_real_, hr_ci_high = NA_real_,
      std_error = NA_real_, statistic = NA_real_, p_value = NA_real_,
      n = nrow(dat), events = sum(dat$os_event == 1),
      model_status = "fit_error", warning_text = conditionMessage(fit)
    ))
  }

  coefficient <- summary(fit)$coefficients["expr", ]
  log_hr <- unname(coefficient[["coef"]])
  std_error <- unname(coefficient[["se(coef)"]])
  ci_low <- log_hr - stats::qnorm(0.975) * std_error
  ci_high <- log_hr + stats::qnorm(0.975) * std_error

  tibble(
    gene_id = gene_id,
    log_hr = log_hr,
    hr = exp(log_hr),
    log_hr_ci_low = ci_low,
    log_hr_ci_high = ci_high,
    hr_ci_low = exp(ci_low),
    hr_ci_high = exp(ci_high),
    std_error = std_error,
    statistic = unname(coefficient[["z"]]),
    p_value = unname(coefficient[["Pr(>|z|)"]]),
    n = fit$n,
    events = fit$nevent,
    model_status = "ok",
    warning_text = if (length(warning_messages) == 0) NA_character_ else
      paste(unique(warning_messages), collapse = " | ")
  )
}

qc_gene_ids <- intersect(deg$gene_id, rownames(vst_mat))
message("All-gene sensitivity: fitting ", length(qc_gene_ids), " QC-filtered genes.")
survival_results <- bind_rows(lapply(qc_gene_ids, fit_gene))

ensembl_ids <- sub("\\..*$", "", qc_gene_ids)
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) |>
  filter(!is.na(SYMBOL)) |>
  arrange(ENSEMBL, SYMBOL) |>
  distinct(ENSEMBL, .keep_all = TRUE)

read_geo_effect <- function(accession) {
  prefix <- tolower(accession)
  out <- read_csv(
    file.path(DIRS$tables, paste0(prefix, "_limma_tumor_vs_normal.csv")),
    show_col_types = FALSE
  ) |>
    transmute(symbol, log2fc = log2FoldChange, p_value = pvalue, fdr = padj) |>
    distinct(symbol, .keep_all = TRUE)
  names(out)[-1] <- paste0(prefix, c("_log2fc", "_p_value", "_fdr"))
  out
}

gse40435 <- read_geo_effect("GSE40435")
gse53757 <- read_geo_effect("GSE53757")
primary_repro <- read_csv(
  file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"),
  show_col_types = FALSE
) |>
  select(tcga_gene_id, primary_reproducible_deg = reproducible_deg)

n_modeled <- sum(survival_results$model_status == "ok")
survival_results <- survival_results |>
  mutate(
    global_fdr = if_else(
      model_status == "ok",
      p.adjust(p_value, method = "BH", n = n_modeled),
      NA_real_
    )
  )

output <- deg |>
  filter(gene_id %in% qc_gene_ids) |>
  transmute(
    gene_id,
    ensembl = sub("\\..*$", "", gene_id),
    base_mean = baseMean,
    tcga_mle_log2fc = log2FoldChange,
    tcga_apeglm_log2fc = log2FoldChange_apeglm,
    tcga_apeglm_lfc_se = lfcSE_apeglm,
    tcga_de_p_value = pvalue,
    tcga_de_fdr = padj,
    primary_tcga_deg_gate = significant
  ) |>
  left_join(symbol_map, by = c("ensembl" = "ENSEMBL")) |>
  rename(symbol = SYMBOL) |>
  left_join(gse40435, by = "symbol") |>
  left_join(gse53757, by = "symbol") |>
  left_join(primary_repro, by = c("gene_id" = "tcga_gene_id")) |>
  mutate(
    primary_reproducible_deg = coalesce(primary_reproducible_deg, FALSE),
    same_direction_tcga_gse40435 = !is.na(gse40435_log2fc) &
      sign(tcga_apeglm_log2fc) == sign(gse40435_log2fc),
    same_direction_tcga_gse53757 = !is.na(gse53757_log2fc) &
      sign(tcga_apeglm_log2fc) == sign(gse53757_log2fc),
    same_direction_all_three = same_direction_tcga_gse40435 & same_direction_tcga_gse53757
  ) |>
  left_join(survival_results, by = "gene_id") |>
  mutate(global_fdr_lt_0_05 = !is.na(global_fdr) & global_fdr < 0.05) |>
  arrange(global_fdr, desc(abs(tcga_apeglm_log2fc)))

summary_output <- tibble(
  metric = c(
    "qc_filtered_genes",
    "successfully_modeled_genes",
    "failed_or_zero_variance_models",
    "global_survival_fdr_lt_0_05",
    "genes_with_both_geo_effects",
    "same_apeglm_direction_all_three_cohorts"
  ),
  value = c(
    length(qc_gene_ids),
    n_modeled,
    length(qc_gene_ids) - n_modeled,
    sum(output$global_fdr_lt_0_05, na.rm = TRUE),
    sum(!is.na(output$gse40435_log2fc) & !is.na(output$gse53757_log2fc)),
    sum(output$same_direction_all_three, na.rm = TRUE)
  )
)

write_csv_atomic(output, FILES$tcga_apeglm_survival)
write_csv_atomic(summary_output, FILES$tcga_apeglm_survival_summary)

message(
  "All-gene apeglm sensitivity complete. Global survival FDR < 0.05: ",
  sum(output$global_fdr_lt_0_05, na.rm = TRUE), " / ", n_modeled
)
