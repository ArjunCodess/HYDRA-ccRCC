source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(survival)
  library(broom)
})

purity_path <- file.path(DIRS$raw, "aran2015_tcga_tumor_purity.xlsx")
if (!file.exists(purity_path)) {
  message("Downloading published TCGA consensus purity estimates.")
  download.file(SOURCE_URLS$aran2015_purity_data, purity_path, mode = "wb", quiet = FALSE)
}

purity_raw <- read_excel(purity_path, skip = 3)
required_purity_columns <- c("Sample ID", "Cancer type", "CPE")
if (!all(required_purity_columns %in% names(purity_raw))) {
  stop(
    "Unexpected Aran et al. purity schema. Missing: ",
    paste(setdiff(required_purity_columns, names(purity_raw)), collapse = ", ")
  )
}

purity <- purity_raw |>
  transmute(
    sample_id = .data[["Sample ID"]],
    cancer_type = .data[["Cancer type"]],
    cpe = suppressWarnings(as.numeric(.data[["CPE"]]))
  ) |>
  filter(cancer_type == "KIRC", !is.na(sample_id), is.finite(cpe)) |>
  mutate(sample_key = substr(sample_id, 1, 16)) |>
  distinct(sample_key, .keep_all = TRUE)

vst_mat <- read_required_rds(FILES$tcga_vst)
coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
candidates <- read_csv(
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  show_col_types = FALSE
) |>
  filter(high_confidence_candidate, tcga_gene_id %in% rownames(vst_mat)) |>
  select(symbol, tcga_gene_id, original_main_log_hr = main_log_hr)

clinical_surv <- clinical |>
  transmute(
    patient_barcode = submitter_id,
    os_time,
    os_event,
    age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
    sex = factor(gender),
    stage = factor(normalize_stage(ajcc_pathologic_stage)),
    grade = factor(case_when(
      normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
      normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
      TRUE ~ NA_character_
    ))
  )

sample_data <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(
    patient_barcode = tcga_patient_barcode(sample_barcode),
    sample_key = substr(sample_barcode, 1, 16)
  ) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  left_join(purity |> select(sample_key, cpe), by = "sample_key") |>
  filter(
    !is.na(os_time),
    os_time > 0,
    !is.na(os_event),
    !is.na(age),
    !is.na(sex),
    !is.na(stage),
    !is.na(grade),
    sample_barcode %in% colnames(vst_mat)
  ) |>
  distinct(patient_barcode, .keep_all = TRUE) |>
  mutate(across(c(sex, stage, grade), droplevels))

fit_candidate <- function(gene_id) {
  dat <- sample_data |>
    mutate(
      expr = as.numeric(scale(vst_mat[gene_id, sample_barcode])),
      purity = as.numeric(scale(cpe))
    ) |>
    select(os_time, os_event, expr, purity, age, sex, stage, grade) |>
    filter(if_all(everything(), ~ !is.na(.x))) |>
    mutate(across(c(sex, stage, grade), droplevels))

  if (nrow(dat) < 100 || sum(dat$os_event) < 25 ||
      any(vapply(dat[c("sex", "stage", "grade")], nlevels, integer(1)) < 2)) {
    return(NULL)
  }

  fit <- tryCatch(
    coxph(
      Surv(os_time, os_event) ~ expr + purity + age + sex + stage + grade,
      data = dat
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  terms <- tidy(fit, conf.int = TRUE)
  gene_term <- terms |> filter(term == "expr")
  purity_term <- terms |> filter(term == "purity")
  if (nrow(gene_term) != 1 || nrow(purity_term) != 1) return(NULL)

  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  gene_ph_p <- if (is.null(zph) || !"expr" %in% rownames(zph$table)) {
    NA_real_
  } else {
    zph$table["expr", "p"]
  }

  tibble(
    gene_log_hr = gene_term$estimate,
    gene_hr = exp(gene_term$estimate),
    gene_hr_ci_low = exp(gene_term$conf.low),
    gene_hr_ci_high = exp(gene_term$conf.high),
    gene_p_value = gene_term$p.value,
    gene_ph_p_value = gene_ph_p,
    purity_log_hr = purity_term$estimate,
    purity_hr = exp(purity_term$estimate),
    purity_p_value = purity_term$p.value,
    n = fit$n,
    events = fit$nevent
  )
}

results <- bind_rows(lapply(seq_len(nrow(candidates)), function(i) {
  out <- fit_candidate(candidates$tcga_gene_id[i])
  if (is.null(out)) return(NULL)
  bind_cols(candidates[i, ], out)
})) |>
  mutate(
    gene_fdr = p.adjust(gene_p_value, method = "BH"),
    same_direction_after_purity = sign(gene_log_hr) == sign(original_main_log_hr),
    absolute_log_hr_attenuation = abs(original_main_log_hr) - abs(gene_log_hr),
    relative_log_hr_attenuation = if_else(
      abs(original_main_log_hr) > 0,
      1 - abs(gene_log_hr) / abs(original_main_log_hr),
      NA_real_
    )
  ) |>
  arrange(gene_fdr)

coverage <- tibble(
  metric = c(
    "published_kirc_purity_samples",
    "tcga_primary_tumors",
    "tcga_primary_tumors_with_cpe",
    "complete_case_purity_survival_samples",
    "complete_case_purity_survival_events",
    "frozen_candidates_tested",
    "same_direction_after_purity",
    "purity_adjusted_fdr_candidates"
  ),
  value = c(
    nrow(purity),
    nrow(sample_data),
    sum(is.finite(sample_data$cpe)),
    if (nrow(results) > 0) unique(results$n)[1] else 0,
    if (nrow(results) > 0) unique(results$events)[1] else 0,
    nrow(results),
    sum(results$same_direction_after_purity, na.rm = TRUE),
    sum(results$gene_fdr < 0.05, na.rm = TRUE)
  )
)

write_csv_atomic(
  results,
  file.path(DIRS$tables, "candidate_direct_tumor_purity_sensitivity.csv")
)
write_csv_atomic(
  coverage,
  file.path(DIRS$tables, "tumor_purity_coverage.csv")
)

message(
  "Direct tumor-purity sensitivity complete. Candidates tested: ",
  nrow(results),
  "; same direction: ",
  sum(results$same_direction_after_purity, na.rm = TRUE),
  "; FDR < 0.05: ",
  sum(results$gene_fdr < 0.05, na.rm = TRUE)
)
