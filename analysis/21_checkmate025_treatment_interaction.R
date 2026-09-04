source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(survival)
  library(tibble)
})

raw_dir <- file.path(DIRS$raw, "checkmate_braun")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
workbook_path <- file.path(raw_dir, "41591_2020_839_MOESM2_ESM.xlsx")

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

download_atomic(
  SOURCE_URLS$checkmate_braun,
  workbook_path,
  identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")
)

workbook_md5 <- unname(tools::md5sum(workbook_path))
if (!identical(tolower(workbook_md5), CHECKMATE_BRAUN_WORKBOOK_MD5)) {
  stop(
    "CheckMate supplementary workbook checksum differs from the pinned value. Expected ",
    CHECKMATE_BRAUN_WORKBOOK_MD5, "; found ", workbook_md5, ".",
    call. = FALSE
  )
}

write_csv_atomic(
  tibble(
    source_file = basename(workbook_path),
    doi = "10.1038/s41591-020-0839-y",
    url = SOURCE_URLS$checkmate_braun,
    bytes = as.numeric(file.info(workbook_path)$size),
    md5 = workbook_md5
  ),
  file.path(DIRS$tables, "checkmate025_source_file.csv")
)

clinical_raw <- read_excel(
  workbook_path,
  sheet = "S1_Clinical_and_Immune_Data",
  skip = 1,
  na = c("", "NA")
)
required_clinical <- c(
  "SUBJID", "Cohort", "Arm", "RNA_ID", "Sex", "Age", "MSKCC",
  "PFS", "PFS_CNSR", "OS", "OS_CNSR",
  "Tumor_Sample_Primary_or_Metastasis"
)
missing_clinical <- setdiff(required_clinical, names(clinical_raw))
if (length(missing_clinical) > 0) {
  stop(
    "CheckMate clinical sheet is missing columns: ",
    paste(missing_clinical, collapse = ", "),
    call. = FALSE
  )
}

clinical <- clinical_raw |>
  transmute(
    patient = as.character(SUBJID),
    cohort = as.character(Cohort),
    arm = toupper(as.character(Arm)),
    rna_id = as.character(RNA_ID),
    sex = factor(as.character(Sex)),
    age = suppressWarnings(as.numeric(Age)),
    mskcc = factor(as.character(MSKCC), levels = c("FAVORABLE", "INTERMEDIATE", "POOR")),
    sample_origin = factor(as.character(Tumor_Sample_Primary_or_Metastasis)),
    pfs_time_months = suppressWarnings(as.numeric(PFS)),
    pfs_event = suppressWarnings(as.integer(PFS_CNSR)),
    os_time_months = suppressWarnings(as.numeric(OS)),
    os_event = suppressWarnings(as.integer(OS_CNSR)),
    treatment_nivolumab = as.integer(arm == "NIVOLUMAB")
  ) |>
  filter(
    cohort == "CM-025",
    arm %in% c("NIVOLUMAB", "EVEROLIMUS"),
    !is.na(rna_id),
    rna_id != "",
    is.finite(pfs_time_months),
    pfs_time_months > 0,
    pfs_event %in% c(0L, 1L),
    is.finite(os_time_months),
    os_time_months > 0,
    os_event %in% c(0L, 1L)
  )

if (nrow(clinical) != 250 || anyDuplicated(clinical$rna_id)) {
  stop("Expected 250 unique RNA-linked CheckMate 025 patients.", call. = FALSE)
}
if (!identical(sort(unique(clinical$arm)), c("EVEROLIMUS", "NIVOLUMAB"))) {
  stop("CheckMate 025 treatment arms were not parsed as expected.", call. = FALSE)
}

candidates <- read_csv(
  file.path(DIRS$tables, "manuscript_candidate_prioritization.csv"),
  show_col_types = FALSE
) |>
  select(symbol, manual_tier, manuscript_role)

expression_raw <- read_excel(
  workbook_path,
  sheet = "S4A_RNA_Expression",
  skip = 1,
  na = c("", "NA"),
  .name_repair = "minimal"
)
if (!"gene_name" %in% names(expression_raw)) {
  stop("CheckMate expression sheet is missing gene_name.", call. = FALSE)
}
missing_samples <- setdiff(clinical$rna_id, names(expression_raw))
if (length(missing_samples) > 0) {
  stop(
    "CheckMate expression matrix is missing clinical RNA IDs: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}

expression_candidates <- expression_raw |>
  filter(gene_name %in% candidates$symbol) |>
  select(gene_name, all_of(clinical$rna_id)) |>
  mutate(across(-gene_name, as.numeric)) |>
  group_by(gene_name) |>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

missing_candidates <- setdiff(candidates$symbol, expression_candidates$gene_name)
if (length(missing_candidates) > 0) {
  stop(
    "CheckMate expression matrix is missing candidates: ",
    paste(missing_candidates, collapse = ", "),
    call. = FALSE
  )
}

term_result <- function(fit, term) {
  if (!term %in% names(coef(fit))) {
    stop("Expected Cox term was not estimable: ", term, call. = FALSE)
  }
  log_hr <- unname(coef(fit)[term])
  standard_error <- sqrt(unname(vcov(fit)[term, term]))
  z_value <- log_hr / standard_error
  ci <- log_hr + qnorm(c(0.025, 0.975)) * standard_error
  ph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph_p <- if (is.null(ph) || !term %in% rownames(ph$table)) {
    NA_real_
  } else {
    unname(ph$table[term, "p"])
  }
  tibble(
    log_hr = log_hr,
    hr = exp(log_hr),
    hr_ci_low = exp(ci[1]),
    hr_ci_high = exp(ci[2]),
    p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
    ph_p_value = ph_p
  )
}

fit_candidate <- function(symbol, endpoint, adjusted) {
  expression_row <- expression_candidates |>
    filter(gene_name == symbol)
  values <- as.numeric(expression_row[1, clinical$rna_id])
  data <- clinical |>
    mutate(expression_z = as.numeric(scale(values)))

  if (adjusted) {
    data <- data |>
      filter(!is.na(age), !is.na(sex), !is.na(mskcc))
  }
  if (any(!is.finite(data$expression_z))) {
    stop("Non-finite standardized expression for ", symbol, ".", call. = FALSE)
  }

  time_column <- paste0(endpoint, "_time_months")
  event_column <- paste0(endpoint, "_event")
  model_data <- data |>
    transmute(
      time = .data[[time_column]],
      event = .data[[event_column]],
      treatment_nivolumab,
      expression_z,
      age,
      sex,
      mskcc
    )
  formula <- if (adjusted) {
    Surv(time, event) ~ treatment_nivolumab * expression_z + age + sex + mskcc
  } else {
    Surv(time, event) ~ treatment_nivolumab * expression_z
  }
  fit <- coxph(formula, data = model_data, ties = "efron")
  coefficient_names <- names(coef(fit))
  interaction_term <- coefficient_names[
    grepl("treatment_nivolumab", coefficient_names, fixed = TRUE) &
      grepl("expression_z", coefficient_names, fixed = TRUE) &
      grepl(":", coefficient_names, fixed = TRUE)
  ]
  if (length(interaction_term) != 1) {
    stop("Candidate-by-treatment interaction term was not uniquely estimable.", call. = FALSE)
  }
  interaction <- term_result(fit, interaction_term)

  arm_effect <- function(arm_value) {
    arm_data <- model_data |>
      filter(treatment_nivolumab == arm_value)
    arm_fit <- coxph(Surv(time, event) ~ expression_z, data = arm_data, ties = "efron")
    term_result(arm_fit, "expression_z")
  }
  nivolumab <- arm_effect(1L)
  everolimus <- arm_effect(0L)

  tibble(
    endpoint = toupper(endpoint),
    model = if_else(adjusted, "age_sex_mskcc_adjusted", "randomized_unadjusted"),
    symbol = symbol,
    n = fit$n,
    events = fit$nevent,
    interaction_log_hr = interaction$log_hr,
    interaction_hr = interaction$hr,
    interaction_hr_ci_low = interaction$hr_ci_low,
    interaction_hr_ci_high = interaction$hr_ci_high,
    interaction_p_value = interaction$p_value,
    interaction_ph_p_value = interaction$ph_p_value,
    nivolumab_log_hr = nivolumab$log_hr,
    nivolumab_p_value = nivolumab$p_value,
    everolimus_log_hr = everolimus$log_hr,
    everolimus_p_value = everolimus$p_value
  )
}

interaction_results <- bind_rows(lapply(c("os", "pfs"), function(endpoint) {
  bind_rows(lapply(c(FALSE, TRUE), function(adjusted) {
    bind_rows(lapply(candidates$symbol, fit_candidate, endpoint = endpoint, adjusted = adjusted))
  }))
})) |>
  group_by(endpoint, model) |>
  mutate(interaction_fdr = p.adjust(interaction_p_value, method = "BH")) |>
  ungroup() |>
  left_join(candidates, by = "symbol") |>
  arrange(endpoint, model, interaction_p_value, symbol)

write_csv_atomic(
  interaction_results,
  file.path(DIRS$tables, "checkmate025_candidate_treatment_interactions.csv")
)

overall_treatment <- bind_rows(lapply(c("os", "pfs"), function(endpoint) {
  time_column <- paste0(endpoint, "_time_months")
  event_column <- paste0(endpoint, "_event")
  data <- clinical |>
    transmute(
      time = .data[[time_column]],
      event = .data[[event_column]],
      treatment_nivolumab
    )
  fit <- coxph(Surv(time, event) ~ treatment_nivolumab, data = data, ties = "efron")
  result <- term_result(fit, "treatment_nivolumab")
  tibble(
    endpoint = toupper(endpoint),
    n = fit$n,
    events = fit$nevent,
    nivolumab_vs_everolimus_hr = result$hr,
    hr_ci_low = result$hr_ci_low,
    hr_ci_high = result$hr_ci_high,
    p_value = result$p_value
  )
}))

study_summary <- tibble(
  metric = c(
    "checkmate025_rna_linked_patients",
    "nivolumab_patients",
    "everolimus_patients",
    "os_events",
    "pfs_events",
    "adjusted_complete_cases",
    "candidates_mapped",
    "os_unadjusted_nominal_interactions",
    "os_unadjusted_fdr_interactions",
    "os_adjusted_nominal_interactions",
    "os_adjusted_fdr_interactions",
    "pfs_unadjusted_nominal_interactions",
    "pfs_unadjusted_fdr_interactions",
    "pfs_adjusted_nominal_interactions",
    "pfs_adjusted_fdr_interactions"
  ),
  value = c(
    nrow(clinical),
    sum(clinical$arm == "NIVOLUMAB"),
    sum(clinical$arm == "EVEROLIMUS"),
    sum(clinical$os_event),
    sum(clinical$pfs_event),
    sum(complete.cases(clinical[, c("age", "sex", "mskcc")])),
    nrow(candidates),
    sum(interaction_results$endpoint == "OS" &
          interaction_results$model == "randomized_unadjusted" &
          interaction_results$interaction_p_value < 0.05),
    sum(interaction_results$endpoint == "OS" &
          interaction_results$model == "randomized_unadjusted" &
          interaction_results$interaction_fdr < 0.05),
    sum(interaction_results$endpoint == "OS" &
          interaction_results$model == "age_sex_mskcc_adjusted" &
          interaction_results$interaction_p_value < 0.05),
    sum(interaction_results$endpoint == "OS" &
          interaction_results$model == "age_sex_mskcc_adjusted" &
          interaction_results$interaction_fdr < 0.05),
    sum(interaction_results$endpoint == "PFS" &
          interaction_results$model == "randomized_unadjusted" &
          interaction_results$interaction_p_value < 0.05),
    sum(interaction_results$endpoint == "PFS" &
          interaction_results$model == "randomized_unadjusted" &
          interaction_results$interaction_fdr < 0.05),
    sum(interaction_results$endpoint == "PFS" &
          interaction_results$model == "age_sex_mskcc_adjusted" &
          interaction_results$interaction_p_value < 0.05),
    sum(interaction_results$endpoint == "PFS" &
          interaction_results$model == "age_sex_mskcc_adjusted" &
          interaction_results$interaction_fdr < 0.05)
  )
)

write_csv_atomic(
  overall_treatment,
  file.path(DIRS$tables, "checkmate025_overall_treatment_effects.csv")
)
write_csv_atomic(
  study_summary,
  file.path(DIRS$tables, "checkmate025_study_summary.csv")
)

message(
  "CheckMate 025 candidate-by-treatment interaction analysis complete. Patients: ",
  nrow(clinical), "; OS events: ", sum(clinical$os_event),
  "; PFS events: ", sum(clinical$pfs_event), "."
)
