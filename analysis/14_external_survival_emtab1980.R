source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readxl)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(survival)
  library(broom)
})

accession <- "E-MTAB-1980"
expression_url <- paste0(
  "https://www.ebi.ac.uk/biostudies/files/",
  accession,
  "/ccRCC_exp_log_quantile_normalized.txt"
)
clinical_url <- paste0(
  "https://static-content.springer.com/esm/",
  "art%3A10.1038%2Fng.2699/MediaObjects/",
  "41588_2013_BFng2699_MOESM35_ESM.xlsx"
)

expression_path <- file.path(DIRS$raw, "arrayexpress", accession, "ccRCC_exp_log_quantile_normalized.txt")
clinical_path <- file.path(DIRS$raw, "arrayexpress", accession, "supplementary_table_1.xlsx")

download_if_needed <- function(url, destination) {
  if (file.exists(destination) && file.info(destination)$size > 0) return(invisible(destination))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(destination, ".download")
  on.exit(unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", quiet = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Downloaded file is empty: ", url, call. = FALSE)
  }
  if (!file.rename(tmp, destination)) {
    stop("Could not move downloaded file into place: ", destination, call. = FALSE)
  }
  invisible(destination)
}

download_if_needed(expression_url, expression_path)
download_if_needed(clinical_url, clinical_path)

clinical_raw <- readxl::read_excel(clinical_path, sheet = 1, skip = 1)
required_clinical_columns <- c(
  "sample ID",
  "Sex",
  "Age",
  "Stage at diagnosis",
  "Fuhrman grade",
  "outcome",
  "observation period (month)",
  "gene expression profile"
)
missing_clinical_columns <- setdiff(required_clinical_columns, names(clinical_raw))
if (length(missing_clinical_columns) > 0) {
  stop(
    "E-MTAB-1980 clinical file is missing columns: ",
    paste(missing_clinical_columns, collapse = ", "),
    call. = FALSE
  )
}

clinical <- clinical_raw |>
  transmute(
    sample = as.character(.data[["sample ID"]]),
    sex = factor(toupper(as.character(.data[["Sex"]]))),
    age = suppressWarnings(as.numeric(.data[["Age"]])),
    pathologic_stage = as.character(.data[["Stage at diagnosis"]]),
    grade = suppressWarnings(as.numeric(.data[["Fuhrman grade"]])),
    molecular_subtype = as.character(.data[["gene expression profile"]]),
    os_event = as.integer(str_to_lower(as.character(.data[["outcome"]])) == "dead"),
    os_time_months = suppressWarnings(as.numeric(.data[["observation period (month)"]])),
    t_high = str_detect(pathologic_stage, "^pT(?:3|4)"),
    node_positive = str_detect(pathologic_stage, "N(?:1|2)"),
    metastatic = str_detect(pathologic_stage, "M1"),
    grade_high = grade >= 3
  ) |>
  filter(
    molecular_subtype %in% c("ccA", "ccB"),
    !is.na(os_time_months),
    os_time_months > 0,
    !is.na(os_event)
  )

if (nrow(clinical) != 101) {
  stop("Expected 101 E-MTAB-1980 expression-profile patients; found ", nrow(clinical), ".", call. = FALSE)
}

expression_raw <- read_tsv(
  expression_path,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

required_expression_columns <- c("REF", "SystematicName")
missing_expression_columns <- setdiff(required_expression_columns, names(expression_raw))
if (length(missing_expression_columns) > 0) {
  stop(
    "E-MTAB-1980 expression file is missing columns: ",
    paste(missing_expression_columns, collapse = ", "),
    call. = FALSE
  )
}

sample_columns <- intersect(clinical$sample, names(expression_raw))
if (length(sample_columns) != nrow(clinical)) {
  missing_samples <- setdiff(clinical$sample, sample_columns)
  stop(
    "E-MTAB-1980 expression matrix is missing clinical samples: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}

refseq_ids <- expression_raw |>
  transmute(refseq = sub("\\.\\d+$", "", as.character(SystematicName))) |>
  filter(!is.na(refseq), refseq != "") |>
  distinct() |>
  pull(refseq)

refseq_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = refseq_ids,
  keytype = "REFSEQ",
  columns = "SYMBOL"
) |>
  filter(!is.na(SYMBOL)) |>
  distinct(REFSEQ, SYMBOL)

expression_gene <- expression_raw |>
  mutate(refseq = sub("\\.\\d+$", "", as.character(SystematicName))) |>
  inner_join(refseq_map, by = c("refseq" = "REFSEQ")) |>
  select(SYMBOL, all_of(sample_columns)) |>
  mutate(across(all_of(sample_columns), as.numeric)) |>
  group_by(SYMBOL) |>
  summarise(across(all_of(sample_columns), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

candidate_priority <- read_csv(
  file.path(DIRS$tables, "manuscript_candidate_prioritization.csv"),
  show_col_types = FALSE
)

fit_external_gene <- function(symbol, tcga_log_hr, manual_tier) {
  row <- expression_gene |> filter(SYMBOL == symbol)
  if (nrow(row) == 0) {
    return(tibble(
      symbol = symbol,
      manual_tier = manual_tier,
      external_present = FALSE,
      external_n = nrow(clinical),
      external_events = sum(clinical$os_event),
      external_log_hr = NA_real_,
      external_hr = NA_real_,
      external_hr_ci_low = NA_real_,
      external_hr_ci_high = NA_real_,
      external_p_value = NA_real_,
      external_ph_p_value = NA_real_,
      external_adjusted_log_hr = NA_real_,
      external_adjusted_p_value = NA_real_,
      external_adjusted_ph_p_value = NA_real_,
      external_same_direction = FALSE,
      external_adjusted_same_direction = FALSE
    ))
  }

  values <- unlist(row[1, sample_columns], use.names = TRUE)
  values <- as.numeric(values)
  names(values) <- sample_columns
  dat <- clinical |>
    mutate(expr = as.numeric(scale(values[sample])))

  unadjusted_fit <- coxph(Surv(os_time_months, os_event) ~ expr, data = dat)
  unadjusted_term <- tidy(unadjusted_fit, conf.int = TRUE) |> filter(term == "expr")
  unadjusted_ph <- tryCatch(cox.zph(unadjusted_fit), error = function(e) NULL)
  unadjusted_ph_p <- if (
    is.null(unadjusted_ph) ||
      !"expr" %in% rownames(unadjusted_ph$table)
  ) {
    NA_real_
  } else {
    unadjusted_ph$table["expr", "p"]
  }

  adjusted_dat <- dat |>
    filter(
      !is.na(age),
      !is.na(t_high),
      !is.na(metastatic)
    )
  adjusted_fit <- tryCatch(
    coxph(
      Surv(os_time_months, os_event) ~ expr + age + t_high + metastatic,
      data = adjusted_dat
    ),
    error = function(e) NULL
  )
  adjusted_term <- if (is.null(adjusted_fit)) {
    NULL
  } else {
    tidy(adjusted_fit, conf.int = TRUE) |> filter(term == "expr")
  }
  adjusted_ph <- if (is.null(adjusted_fit)) {
    NULL
  } else {
    tryCatch(cox.zph(adjusted_fit), error = function(e) NULL)
  }
  adjusted_ph_p <- if (
    is.null(adjusted_ph) ||
      !"expr" %in% rownames(adjusted_ph$table)
  ) {
    NA_real_
  } else {
    adjusted_ph$table["expr", "p"]
  }

  tibble(
    symbol = symbol,
    manual_tier = manual_tier,
    external_present = TRUE,
    external_n = unadjusted_fit$n,
    external_events = unadjusted_fit$nevent,
    external_log_hr = unadjusted_term$estimate,
    external_hr = exp(unadjusted_term$estimate),
    external_hr_ci_low = exp(unadjusted_term$conf.low),
    external_hr_ci_high = exp(unadjusted_term$conf.high),
    external_p_value = unadjusted_term$p.value,
    external_ph_p_value = unadjusted_ph_p,
    external_adjusted_log_hr = if (is.null(adjusted_term)) NA_real_ else adjusted_term$estimate,
    external_adjusted_p_value = if (is.null(adjusted_term)) NA_real_ else adjusted_term$p.value,
    external_adjusted_ph_p_value = adjusted_ph_p,
    external_same_direction = sign(unadjusted_term$estimate) == sign(tcga_log_hr),
    external_adjusted_same_direction = if (is.null(adjusted_term)) {
      FALSE
    } else {
      sign(adjusted_term$estimate) == sign(tcga_log_hr)
    }
  )
}

external_results <- bind_rows(lapply(seq_len(nrow(candidate_priority)), function(i) {
  fit_external_gene(
    symbol = candidate_priority$symbol[i],
    tcga_log_hr = candidate_priority$main_log_hr[i],
    manual_tier = candidate_priority$manual_tier[i]
  )
})) |>
  mutate(
    external_fdr = p.adjust(external_p_value, method = "BH"),
    external_adjusted_fdr = p.adjust(external_adjusted_p_value, method = "BH"),
    external_directional_nominal_support = external_present &
      external_same_direction &
      !is.na(external_p_value) &
      external_p_value < 0.05,
    external_adjusted_directional_nominal_support = external_present &
      external_adjusted_same_direction &
      !is.na(external_adjusted_p_value) &
      external_adjusted_p_value < 0.05,
    external_strict_support = external_present &
      external_same_direction &
      !is.na(external_fdr) &
      external_fdr < 0.05 &
      !is.na(external_ph_p_value) &
      external_ph_p_value >= THRESHOLDS$ph_min_p &
      external_adjusted_same_direction &
      !is.na(external_adjusted_fdr) &
      external_adjusted_fdr < 0.05 &
      !is.na(external_adjusted_ph_p_value) &
      external_adjusted_ph_p_value >= THRESHOLDS$ph_min_p,
    external_interpretation = case_when(
      !external_present ~ "not mapped in E-MTAB-1980",
      external_strict_support ~
        "same-direction FDR support with limited adjustment and proportional-hazards support",
      external_same_direction &
        !is.na(external_fdr) &
        external_fdr < 0.05 &
        !is.na(external_ph_p_value) &
        external_ph_p_value < THRESHOLDS$ph_min_p ~
        "same-direction association with proportional-hazards violation",
      external_directional_nominal_support ~ "same-direction nominal external support",
      external_same_direction ~ "same direction but not nominally significant",
      !external_same_direction & !is.na(external_p_value) & external_p_value < 0.05 ~
        "opposite-direction nominal external signal",
      TRUE ~ "not externally supported in this cohort"
    )
  ) |>
  left_join(
    candidate_priority |>
      select(
        symbol,
        main_log_hr,
        main_hr,
        main_fdr,
        manuscript_role,
        survival_direction,
        tumor_direction
      ),
    by = "symbol"
  ) |>
  arrange(
    factor(
      manual_tier,
      levels = c(
        "lead",
        "supporting",
        "supporting risk",
        "interpret cautiously",
        "composition flag",
        "do not highlight"
      )
    ),
    external_p_value
  )

write_csv_atomic(
  external_results,
  file.path(DIRS$tables, "external_survival_emtab1980.csv")
)

write_csv_atomic(
  tibble(
    metric = c(
      "emtab1980_samples",
      "emtab1980_events",
      "tested_candidates",
      "platform_present_candidates",
      "same_direction_candidates",
      "same_direction_nominal_candidates",
      "same_direction_fdr_candidates",
      "same_direction_adjusted_nominal_candidates",
      "same_direction_adjusted_fdr_candidates",
      "proportional_hazards_pass_candidates",
      "adjusted_proportional_hazards_pass_candidates",
      "strict_external_support_candidates",
      "opposite_direction_nominal_candidates"
    ),
    value = c(
      nrow(clinical),
      sum(clinical$os_event),
      nrow(external_results),
      sum(external_results$external_present),
      sum(external_results$external_same_direction),
      sum(external_results$external_directional_nominal_support),
      sum(
        external_results$external_same_direction &
          !is.na(external_results$external_fdr) &
          external_results$external_fdr < 0.05
      ),
      sum(external_results$external_adjusted_directional_nominal_support),
      sum(
        external_results$external_adjusted_same_direction &
          !is.na(external_results$external_adjusted_fdr) &
          external_results$external_adjusted_fdr < 0.05
      ),
      sum(
        external_results$external_present &
          !is.na(external_results$external_ph_p_value) &
          external_results$external_ph_p_value >= THRESHOLDS$ph_min_p
      ),
      sum(
        external_results$external_present &
          !is.na(external_results$external_adjusted_ph_p_value) &
          external_results$external_adjusted_ph_p_value >= THRESHOLDS$ph_min_p
      ),
      sum(external_results$external_strict_support),
      sum(
        !external_results$external_same_direction &
          !is.na(external_results$external_p_value) &
          external_results$external_p_value < 0.05
      )
    )
  ),
  file.path(DIRS$tables, "external_survival_emtab1980_summary.csv")
)

message(
  "E-MTAB-1980 frozen-candidate external survival test complete. Samples: ",
  nrow(clinical),
  "; events: ",
  sum(clinical$os_event),
  "."
)
