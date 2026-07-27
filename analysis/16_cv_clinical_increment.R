source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

vst_mat <- read_required_rds(FILES$tcga_vst)
coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
candidates <- read_csv(file.path(DIRS$tables, "candidate_gene_evidence_table.csv"), show_col_types = FALSE) |>
  filter(high_confidence_candidate, tcga_gene_id %in% rownames(vst_mat)) |>
  select(symbol, tcga_gene_id, main_log_hr)

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

dat <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  filter(!is.na(os_time), os_time > 0, !is.na(os_event), !is.na(age),
         !is.na(sex), !is.na(stage), !is.na(grade), sample_barcode %in% colnames(vst_mat)) |>
  distinct(patient_barcode, .keep_all = TRUE) |>
  mutate(across(c(sex, stage, grade), droplevels))

make_folds <- function(event, k) {
  folds <- integer(length(event))
  for (value in sort(unique(event))) {
    idx <- which(event == value)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

c_index <- function(time, event, score) {
  unname(concordance(Surv(time, event) ~ score, reverse = TRUE)$concordance)
}

set.seed(RESAMPLING$seed + 1L)
repeat_rows <- vector("list", RESAMPLING$cv_repeats * nrow(candidates))
row_id <- 1L

for (repeat_id in seq_len(RESAMPLING$cv_repeats)) {
  message("Clinical-increment CV repeat ", repeat_id, "/", RESAMPLING$cv_repeats)
  fold_id <- make_folds(dat$os_event, RESAMPLING$cv_folds)

  for (i in seq_len(nrow(candidates))) {
    gene <- candidates[i, ]
    gene_dat <- dat |>
      mutate(expr = as.numeric(scale(vst_mat[gene$tcga_gene_id, sample_barcode])))
    clinical_lp <- rep(NA_real_, nrow(gene_dat))
    gene_lp <- rep(NA_real_, nrow(gene_dat))

    for (fold in seq_len(RESAMPLING$cv_folds)) {
      train <- gene_dat[fold_id != fold, ]
      test <- gene_dat[fold_id == fold, ]
      clinical_fit <- tryCatch(
        coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train, x = TRUE),
        error = function(e) NULL
      )
      gene_fit <- tryCatch(
        coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade, data = train, x = TRUE),
        error = function(e) NULL
      )
      if (is.null(clinical_fit) || is.null(gene_fit)) next
      clinical_lp[fold_id == fold] <- tryCatch(
        predict(clinical_fit, newdata = test, type = "lp", reference = "zero"),
        error = function(e) rep(NA_real_, nrow(test))
      )
      gene_lp[fold_id == fold] <- tryCatch(
        predict(gene_fit, newdata = test, type = "lp", reference = "zero"),
        error = function(e) rep(NA_real_, nrow(test))
      )
    }

    complete <- is.finite(clinical_lp) & is.finite(gene_lp)
    if (sum(complete) < 100) next
    clinical_c <- c_index(gene_dat$os_time[complete], gene_dat$os_event[complete], clinical_lp[complete])
    gene_c <- c_index(gene_dat$os_time[complete], gene_dat$os_event[complete], gene_lp[complete])
    repeat_rows[[row_id]] <- tibble(
      repeat_id,
      symbol = gene$symbol,
      n_tested = sum(complete),
      clinical_c_index = clinical_c,
      clinical_gene_c_index = gene_c,
      delta_c_index = gene_c - clinical_c
    )
    row_id <- row_id + 1L
  }
}

repeat_table <- bind_rows(repeat_rows)
summary_table <- repeat_table |>
  group_by(symbol) |>
  summarise(
    repeats = n(),
    mean_clinical_c_index = mean(clinical_c_index),
    mean_clinical_gene_c_index = mean(clinical_gene_c_index),
    mean_delta_c_index = mean(delta_c_index),
    median_delta_c_index = median(delta_c_index),
    delta_ci_low = quantile(delta_c_index, 0.025),
    delta_ci_high = quantile(delta_c_index, 0.975),
    positive_delta_frequency = mean(delta_c_index > 0),
    .groups = "drop"
  ) |>
  left_join(candidates |> select(symbol, main_log_hr), by = "symbol") |>
  arrange(desc(mean_delta_c_index))

write_csv_atomic(repeat_table, file.path(DIRS$tables, "candidate_cv_clinical_increment_repeats.csv"))
write_csv_atomic(summary_table, file.path(DIRS$tables, "candidate_cv_clinical_increment.csv"))

message("Clinical-increment CV complete. Candidates assessed: ", nrow(summary_table))
