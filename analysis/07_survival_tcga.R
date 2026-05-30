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
    stage_clean = normalize_stage(ajcc_pathologic_stage),
    grade_clean = normalize_grade(tumor_grade)
  )

tumor_samples <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  filter(!is.na(os_time), os_time > 0, !is.na(os_event))

repro_path <- file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv")
if (file.exists(repro_path)) {
  candidate_genes <- read_csv(repro_path, show_col_types = FALSE) |>
    filter(reproducible_deg) |>
    pull(tcga_gene_id)
} else {
  candidate_genes <- deg |>
    filter(!is.na(padj), padj < THRESHOLDS$deg_fdr, abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc) |>
    pull(gene_id)
}

candidate_genes <- intersect(candidate_genes, rownames(vst_mat))
sample_ids <- intersect(tumor_samples$sample_barcode, colnames(vst_mat))
tumor_samples <- tumor_samples |> filter(sample_barcode %in% sample_ids)

fit_gene_model <- function(gene_id, model_type, covars) {
  expr <- as.numeric(scale(vst_mat[gene_id, tumor_samples$sample_barcode]))
  dat_full <- tumor_samples |>
    transmute(
      os_time = os_time,
      os_event = os_event,
      expr = expr,
      age = age_years,
      sex = factor(sex),
      stage = factor(stage_clean),
      grade = factor(case_when(
        grade_clean %in% c("G1", "G2") ~ "Low grade",
        grade_clean %in% c("G3", "G4") ~ "High grade",
        TRUE ~ NA_character_
      ))
    )

  dat <- dat_full |>
    select(os_time, os_event, all_of(covars)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  dat <- dat |>
    mutate(across(any_of(c("sex", "stage", "grade")), droplevels))

  factor_covars <- intersect(c("sex", "stage", "grade"), covars)
  if (any(vapply(dat[factor_covars], nlevels, integer(1)) < 2)) return(NULL)
  if (nrow(dat) < 100 || sum(dat$os_event == 1) < 25) return(NULL)

  form <- as.formula(paste("Surv(os_time, os_event) ~", paste(covars, collapse = " + ")))
  warning_messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxph(form, data = dat),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  ph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph_p <- if (is.null(ph) || !"expr" %in% rownames(ph$table)) NA_real_ else ph$table["expr", "p"]

  tidy(fit, conf.int = TRUE) |>
    filter(term == "expr") |>
    transmute(
      gene_id = gene_id,
      model_type = model_type,
      log_hr = estimate,
      hr = exp(estimate),
      log_hr_ci_low = conf.low,
      log_hr_ci_high = conf.high,
      hr_ci_low = exp(conf.low),
      hr_ci_high = exp(conf.high),
      std_error = std.error,
      statistic = statistic,
      p_value = p.value,
      ph_p_value = ph_p,
      n = fit$n,
      events = fit$nevent,
      covariates = paste(covars, collapse = ";"),
      warning_count = length(warning_messages),
      warning_text = paste(unique(warning_messages), collapse = " | ")
    )
}

model_specs <- list(
  stage_grade_complete = c("expr", "age", "sex", "stage", "grade"),
  stage_complete = c("expr", "age", "sex", "stage"),
  grade_complete = c("expr", "age", "sex", "grade")
)

surv_results <- bind_rows(lapply(names(model_specs), function(model_type) {
  bind_rows(lapply(candidate_genes, fit_gene_model, model_type = model_type, covars = model_specs[[model_type]]))
})) |>
  group_by(model_type) |>
  mutate(fdr = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(fdr)

write_csv_atomic(surv_results, FILES$tcga_survival)

message("Cox modeling complete. Model rows: ", nrow(surv_results))
message("Main model genes FDR < ", THRESHOLDS$strict_survival_fdr, ": ",
        sum(surv_results$model_type == "stage_grade_complete" &
              surv_results$fdr < THRESHOLDS$strict_survival_fdr, na.rm = TRUE))
