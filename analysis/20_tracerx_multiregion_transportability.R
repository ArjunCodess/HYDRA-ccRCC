source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(survival)
})

raw_dir <- file.path(DIRS$raw, "tracerx_renal")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- tibble(
  source_file = c("tumour_tpm.rds", "tx_annotation.tsv", "TRACERx_s1_1_clinical.txt"),
  url = c(
    SOURCE_URLS$tracerx_tumour_tpm,
    SOURCE_URLS$tracerx_annotation,
    SOURCE_URLS$tracerx_clinical
  ),
  path = file.path(raw_dir, source_file)
)

download_atomic <- function(url, destination, force = FALSE) {
  if (!force && file.exists(destination) && file.info(destination)$size > 0) {
    return(invisible(destination))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(destination, ".download")
  on.exit(unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", quiet = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Downloaded file is empty: ", url, call. = FALSE)
  }
  if (file.exists(destination) && unlink(destination) != 0) {
    stop("Could not replace downloaded file: ", destination, call. = FALSE)
  }
  if (!file.rename(tmp, destination)) {
    stop("Could not move downloaded file into place: ", destination, call. = FALSE)
  }
  invisible(destination)
}

force_download <- identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")
for (i in seq_len(nrow(source_files))) {
  download_atomic(source_files$url[i], source_files$path[i], force_download)
}

source_inventory <- source_files |>
  mutate(
    commit = TRACERX_DATA_COMMIT,
    bytes = as.numeric(file.info(path)$size),
    md5 = unname(tools::md5sum(path))
  ) |>
  select(source_file, commit, url, bytes, md5)

write_csv_atomic(
  source_inventory,
  file.path(DIRS$tables, "tracerx_multiregion_source_files.csv")
)

expression <- readRDS(file.path(raw_dir, "tumour_tpm.rds"))
if (!is.matrix(expression) || is.null(rownames(expression)) || is.null(colnames(expression))) {
  stop("TRACERx tumor TPM data must be a named gene-by-sample matrix.")
}

annotation <- read_tsv(
  file.path(raw_dir, "tx_annotation.tsv"),
  show_col_types = FALSE,
  progress = FALSE
)
required_annotation <- c("Patient", "sample", "type_collapsed")
if (!all(required_annotation %in% names(annotation))) {
  stop("TRACERx annotation is missing required columns.")
}

clinical_raw <- read_tsv(
  file.path(raw_dir, "TRACERx_s1_1_clinical.txt"),
  skip = 1,
  na = c("", "-", "N/A"),
  show_col_types = FALSE,
  progress = FALSE
)
required_clinical <- c("Subject", "Total follow up (months)", "Outcome")
if (!all(required_clinical %in% names(clinical_raw))) {
  stop("TRACERx clinical table is missing required survival columns.")
}

clinical <- clinical_raw |>
  transmute(
    patient = as.character(Subject),
    os_time_months = suppressWarnings(as.numeric(`Total follow up (months)`)),
    os_event = as.integer(tolower(as.character(Outcome)) == "death")
  ) |>
  filter(is.finite(os_time_months), os_time_months > 0, !is.na(os_event))

if (anyDuplicated(clinical$patient)) {
  stop("TRACERx clinical table contains duplicate patient identifiers.")
}

sample_data <- annotation |>
  transmute(
    patient = as.character(Patient),
    sample = as.character(sample),
    type_collapsed = as.character(type_collapsed)
  ) |>
  filter(type_collapsed == "PRIMARY", sample %in% colnames(expression)) |>
  inner_join(clinical, by = "patient") |>
  distinct(patient, sample, .keep_all = TRUE)

if (anyDuplicated(sample_data$sample)) {
  stop("TRACERx primary-region annotation contains duplicate sample identifiers.")
}
if (n_distinct(sample_data$patient) < RESAMPLING$tracerx_small_cohort_size) {
  stop("TRACERx has too few survival-linked primary tumors for the size-matched analysis.")
}

candidates <- read_csv(
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  show_col_types = FALSE
) |>
  filter(high_confidence_candidate) |>
  select(symbol, tcga_log_hr = main_log_hr)

missing_candidates <- setdiff(candidates$symbol, rownames(expression))
if (length(missing_candidates) > 0) {
  stop(
    "TRACERx expression matrix is missing candidates: ",
    paste(missing_candidates, collapse = ", ")
  )
}

candidate_expression <- log2(expression[candidates$symbol, sample_data$sample, drop = FALSE] + 1)
expression_long <- as.data.frame(t(candidate_expression), check.names = FALSE) |>
  rownames_to_column("sample") |>
  pivot_longer(-sample, names_to = "symbol", values_to = "expression") |>
  inner_join(sample_data, by = "sample")

if (any(!is.finite(expression_long$expression))) {
  stop("TRACERx candidate expression contains non-finite values.")
}

patient_expression <- expression_long |>
  group_by(symbol, patient, os_time_months, os_event) |>
  summarise(
    patient_median_expression = median(expression),
    region_expression_min = min(expression),
    region_expression_max = max(expression),
    region_expression_range = region_expression_max - region_expression_min,
    n_regions = n(),
    .groups = "drop"
  )

thresholds <- patient_expression |>
  group_by(symbol) |>
  summarise(
    readout_threshold = median(patient_median_expression),
    .groups = "drop"
  )

patient_discordance <- expression_long |>
  left_join(thresholds, by = "symbol") |>
  group_by(symbol, patient) |>
  summarise(
    n_regions = n(),
    has_low_readout = any(expression < readout_threshold),
    has_high_readout = any(expression >= readout_threshold),
    discordant_readout = has_low_readout & has_high_readout,
    .groups = "drop"
  )

discordance_summary <- patient_expression |>
  left_join(patient_discordance, by = c("symbol", "patient", "n_regions")) |>
  left_join(thresholds, by = "symbol") |>
  group_by(symbol, readout_threshold) |>
  summarise(
    patients_total = n(),
    multiregion_patients = sum(n_regions >= 2),
    primary_regions = sum(n_regions),
    discordant_multiregion_patients = sum(discordant_readout & n_regions >= 2),
    discordant_multiregion_percent = 100 *
      discordant_multiregion_patients / multiregion_patients,
    median_within_patient_log2_tpm_range = median(region_expression_range[n_regions >= 2]),
    max_within_patient_log2_tpm_range = max(region_expression_range[n_regions >= 2]),
    .groups = "drop"
  ) |>
  left_join(candidates, by = "symbol") |>
  arrange(desc(discordant_multiregion_percent), symbol)

write_csv_atomic(
  patient_discordance,
  file.path(DIRS$tables, "tracerx_candidate_patient_region_discordance.csv")
)
write_csv_atomic(
  discordance_summary,
  file.path(DIRS$tables, "tracerx_candidate_multiregion_summary.csv")
)

fit_gene <- function(values, time, event) {
  expression_sd <- sd(values)
  if (!is.finite(expression_sd) || expression_sd == 0) {
    return(c(log_hr = NA_real_, p_value = NA_real_))
  }
  data <- tibble(
    os_time_months = time,
    os_event = event,
    expression_z = as.numeric(scale(values))
  )
  fit <- tryCatch(
    coxph(Surv(os_time_months, os_event) ~ expression_z, data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(log_hr = NA_real_, p_value = NA_real_))
  coefficient <- unname(coef(fit)[[1]])
  p_value <- unname(summary(fit)$coefficients[1, "Pr(>|z|)"])
  if (!is.finite(coefficient) || !is.finite(p_value)) {
    return(c(log_hr = NA_real_, p_value = NA_real_))
  }
  c(log_hr = coefficient, p_value = p_value)
}

patient_median_matrix <- patient_expression |>
  select(patient, symbol, patient_median_expression) |>
  pivot_wider(names_from = symbol, values_from = patient_median_expression) |>
  inner_join(clinical, by = "patient") |>
  arrange(patient)

reference_rows <- lapply(seq_len(nrow(candidates)), function(i) {
  symbol <- candidates$symbol[i]
  estimate <- fit_gene(
    patient_median_matrix[[symbol]],
    patient_median_matrix$os_time_months,
    patient_median_matrix$os_event
  )
  tibble(
    symbol = symbol,
    patient_median_log_hr = estimate[["log_hr"]],
    patient_median_hr = exp(estimate[["log_hr"]]),
    patient_median_p_value = estimate[["p_value"]]
  )
})
reference_results <- bind_rows(reference_rows)

patient_status <- sample_data |>
  distinct(patient, os_event)
event_patients <- patient_status$patient[patient_status$os_event == 1]
censored_patients <- patient_status$patient[patient_status$os_event == 0]
small_event_count <- round(
  RESAMPLING$tracerx_small_cohort_size *
    length(event_patients) / nrow(patient_status)
)
small_censored_count <- RESAMPLING$tracerx_small_cohort_size - small_event_count
if (small_event_count < 1 || small_event_count > length(event_patients) ||
    small_censored_count > length(censored_patients)) {
  stop("TRACERx event-stratified small-cohort sample cannot be constructed.")
}

set.seed(RESAMPLING$seed + 20L)
repeat_rows <- vector("list", RESAMPLING$tracerx_region_repeats * 2L)
row_index <- 1L

for (repeat_id in seq_len(RESAMPLING$tracerx_region_repeats)) {
  if (repeat_id %% 100 == 0) {
    message(
      "TRACERx one-region resampling repeat ", repeat_id, "/",
      RESAMPLING$tracerx_region_repeats
    )
  }

  small_patients <- c(
    sample(event_patients, small_event_count, replace = FALSE),
    sample(censored_patients, small_censored_count, replace = FALSE)
  )
  scenario_patients <- list(
    full_cohort = patient_status$patient,
    size_matched_39 = small_patients
  )

  for (scenario in names(scenario_patients)) {
    selected_regions <- sample_data |>
      filter(patient %in% scenario_patients[[scenario]]) |>
      group_by(patient) |>
      slice_sample(n = 1) |>
      ungroup() |>
      arrange(patient)

    selected_expression <- candidate_expression[
      candidates$symbol,
      selected_regions$sample,
      drop = FALSE
    ]
    estimates <- lapply(seq_len(nrow(candidates)), function(i) {
      estimate <- fit_gene(
        as.numeric(selected_expression[i, ]),
        selected_regions$os_time_months,
        selected_regions$os_event
      )
      tibble(
        symbol = candidates$symbol[i],
        tcga_log_hr = candidates$tcga_log_hr[i],
        log_hr = estimate[["log_hr"]],
        hr = exp(estimate[["log_hr"]]),
        p_value = estimate[["p_value"]]
      )
    })

    repeat_rows[[row_index]] <- bind_rows(estimates) |>
      mutate(
        scenario = scenario,
        repeat_id = repeat_id,
        n = nrow(selected_regions),
        events = sum(selected_regions$os_event),
        same_tcga_direction = sign(log_hr) == sign(tcga_log_hr)
      ) |>
      select(
        scenario, repeat_id, n, events, symbol, tcga_log_hr,
        log_hr, hr, p_value, same_tcga_direction
      )
    row_index <- row_index + 1L
  }
}

repeat_results <- bind_rows(repeat_rows)
resampling_summary <- repeat_results |>
  group_by(scenario, symbol, tcga_log_hr) |>
  summarise(
    successful_repeats = sum(is.finite(log_hr)),
    n = first(n),
    events = first(events),
    median_log_hr = median(log_hr, na.rm = TRUE),
    log_hr_ci_low = quantile(log_hr, 0.025, na.rm = TRUE),
    log_hr_ci_high = quantile(log_hr, 0.975, na.rm = TRUE),
    log_hr_sd = sd(log_hr, na.rm = TRUE),
    same_tcga_direction_fraction = mean(same_tcga_direction, na.rm = TRUE),
    nominal_p_lt_0_05_fraction = mean(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(reference_results, by = "symbol") |>
  arrange(scenario, desc(same_tcga_direction_fraction), symbol)

study_summary <- tibble(
  metric = c(
    "tracerx_primary_regions",
    "tracerx_survival_linked_patients",
    "tracerx_survival_events",
    "tracerx_multiregion_patients",
    "candidates_mapped",
    "region_resampling_repeats",
    "size_matched_patients",
    "size_matched_events"
  ),
  value = c(
    nrow(sample_data),
    n_distinct(sample_data$patient),
    sum(patient_status$os_event),
    sum(table(sample_data$patient) >= 2),
    nrow(candidates),
    RESAMPLING$tracerx_region_repeats,
    RESAMPLING$tracerx_small_cohort_size,
    small_event_count
  )
)

write_csv_atomic(
  repeat_results,
  file.path(DIRS$tables, "tracerx_one_region_cox_repeats.csv")
)
write_csv_atomic(
  resampling_summary,
  file.path(DIRS$tables, "tracerx_one_region_cox_summary.csv")
)
write_csv_atomic(
  study_summary,
  file.path(DIRS$tables, "tracerx_multiregion_study_summary.csv")
)

message(
  "TRACERx multiregion transportability analysis complete. Primary regions: ",
  nrow(sample_data), "; patients: ", n_distinct(sample_data$patient),
  "; events: ", sum(patient_status$os_event), "."
)
