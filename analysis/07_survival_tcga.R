source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(broom)
})

vst_mat <- read_required_rds(FILES$tcga_vst)
coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
deg <- read_csv(FILES$tcga_deg, show_col_types = FALSE)

clinical_surv <- clinical |>
  transmute(
    patient_barcode = submitter_id,
    os_time = os_time,
    os_event = os_event,
    age_years = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
    sex = gender,
    stage_clean = stage_clean,
    grade_clean = grade_clean
  )

tumor_samples <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  filter(!is.na(os_time), os_time > 0, !is.na(os_event))

candidate_genes <- deg |>
  filter(!is.na(padj), padj < THRESHOLDS$deg_fdr, abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc) |>
  pull(gene_id)

candidate_genes <- intersect(candidate_genes, rownames(vst_mat))
sample_ids <- intersect(tumor_samples$sample_barcode, colnames(vst_mat))
tumor_samples <- tumor_samples |> filter(sample_barcode %in% sample_ids)

fit_gene <- function(gene_id) {
  expr <- as.numeric(scale(vst_mat[gene_id, tumor_samples$sample_barcode]))
  dat_full <- tumor_samples |>
    transmute(
      os_time = os_time,
      os_event = os_event,
      expr = expr,
      age = age_years,
      sex = factor(sex),
      stage = factor(stage_clean),
      grade = factor(grade_clean)
    )

  candidate_covars <- c("expr", "age", "sex", "stage", "grade")
  usable <- candidate_covars[vapply(dat_full[candidate_covars], function(x) {
    non_missing <- x[!is.na(x)]
    length(non_missing) > 0 && length(unique(non_missing)) >= 2
  }, logical(1))]

  if (!"expr" %in% usable) return(NULL)

  dat <- dat_full |>
    select(os_time, os_event, all_of(usable)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  if (nrow(dat) < 50 || length(unique(dat$os_event)) < 2) return(NULL)

  form <- as.formula(paste("Surv(os_time, os_event) ~", paste(usable, collapse = " + ")))
  fit <- tryCatch(coxph(form, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  ph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph_p <- if (is.null(ph) || !"expr" %in% rownames(ph$table)) NA_real_ else ph$table["expr", "p"]

  tidy(fit) |>
    filter(term == "expr") |>
    transmute(
      gene_id = gene_id,
      log_hr = estimate,
      hr = exp(estimate),
      std_error = std.error,
      statistic = statistic,
      p_value = p.value,
      ph_p_value = ph_p,
      n = fit$n,
      events = fit$nevent,
      covariates = paste(usable, collapse = ";")
    )
}

surv_results <- bind_rows(lapply(candidate_genes, fit_gene)) |>
  mutate(fdr = p.adjust(p_value, method = "BH")) |>
  arrange(fdr)

write_csv_atomic(surv_results, FILES$tcga_survival)

message("Cox modeling complete. Tested genes: ", nrow(surv_results))
message("FDR < ", THRESHOLDS$survival_fdr, ": ", sum(surv_results$fdr < THRESHOLDS$survival_fdr, na.rm = TRUE))
