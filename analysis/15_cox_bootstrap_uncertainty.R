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
candidates <- read_csv(
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  show_col_types = FALSE
)

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
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
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

candidate_table <- candidates |>
  filter(high_confidence_candidate, tcga_gene_id %in% rownames(vst_mat)) |>
  select(symbol, tcga_gene_id, full_fit_log_hr = main_log_hr)

draw_stratified_bootstrap <- function(event) {
  unlist(lapply(sort(unique(event)), function(value) {
    indices <- which(event == value)
    sample(indices, length(indices), replace = TRUE)
  }), use.names = FALSE)
}

fit_bootstrap_gene <- function(gene_id, indices) {
  raw_expr <- as.numeric(vst_mat[gene_id, sample_data$sample_barcode[indices]])
  expr_sd <- sd(raw_expr)
  if (!is.finite(expr_sd) || expr_sd <= 0) return(NULL)

  dat <- sample_data[indices, ] |>
    mutate(expr = (raw_expr - mean(raw_expr)) / expr_sd) |>
    mutate(across(c(sex, stage, grade), droplevels))

  fit <- tryCatch(
    suppressWarnings(coxph(
      Surv(os_time, os_event) ~ expr + age + sex + stage + grade,
      data = dat
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  term <- tidy(fit) |> filter(term == "expr")
  if (nrow(term) != 1 || !is.finite(term$estimate)) return(NULL)

  tibble(
    bootstrap_log_hr = term$estimate,
    bootstrap_hr = exp(term$estimate),
    n = fit$n,
    events = fit$nevent
  )
}

set.seed(RESAMPLING$seed + 1L)
repeat_rows <- vector("list", RESAMPLING$coefficient_bootstrap_repeats)

for (repeat_id in seq_len(RESAMPLING$coefficient_bootstrap_repeats)) {
  indices <- draw_stratified_bootstrap(sample_data$os_event)
  if (repeat_id %% 50L == 0L || repeat_id == 1L) {
    message(
      "Cox coefficient bootstrap repeat ",
      repeat_id,
      "/",
      RESAMPLING$coefficient_bootstrap_repeats
    )
  }

  repeat_rows[[repeat_id]] <- bind_rows(lapply(seq_len(nrow(candidate_table)), function(i) {
    result <- fit_bootstrap_gene(candidate_table$tcga_gene_id[i], indices)
    if (is.null(result)) return(NULL)
    mutate(
      result,
      repeat_id = repeat_id,
      symbol = candidate_table$symbol[i],
      tcga_gene_id = candidate_table$tcga_gene_id[i],
      .before = 1
    )
  }))
}

bootstrap_repeats <- bind_rows(repeat_rows)
bootstrap_summary <- bootstrap_repeats |>
  group_by(symbol, tcga_gene_id) |>
  summarise(
    successful_repeats = n(),
    bootstrap_mean_log_hr = mean(bootstrap_log_hr),
    bootstrap_median_log_hr = median(bootstrap_log_hr),
    bootstrap_se = sd(bootstrap_log_hr),
    bootstrap_ci_low = quantile(bootstrap_log_hr, 0.025),
    bootstrap_ci_high = quantile(bootstrap_log_hr, 0.975),
    .groups = "drop"
  ) |>
  left_join(candidate_table, by = c("symbol", "tcga_gene_id")) |>
  mutate(
    bootstrap_bias = bootstrap_mean_log_hr - full_fit_log_hr,
    direction_agreement = sign(bootstrap_median_log_hr) == sign(full_fit_log_hr),
    ci_excludes_zero = bootstrap_ci_low > 0 | bootstrap_ci_high < 0
  ) |>
  arrange(desc(ci_excludes_zero), desc(abs(full_fit_log_hr)))

write_csv_atomic(
  bootstrap_repeats,
  file.path(DIRS$tables, "candidate_cox_bootstrap_repeats.csv")
)
write_csv_atomic(
  bootstrap_summary,
  file.path(DIRS$tables, "candidate_cox_bootstrap_summary.csv")
)

message(
  "Cox coefficient bootstrap complete. Candidates: ",
  nrow(bootstrap_summary),
  "; requested repeats: ",
  RESAMPLING$coefficient_bootstrap_repeats
)
