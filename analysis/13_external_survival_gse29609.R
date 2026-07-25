source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(survival)
  library(broom)
})

collapse_gse29609_to_gene <- function(expr, feature_data) {
  symbol_col <- names(feature_data)[str_detect(tolower(names(feature_data)), "gene symbol")]
  if (length(symbol_col) == 0) {
    stop("No gene-symbol column found for GSE29609 platform annotation.", call. = FALSE)
  }

  symbols <- as.character(feature_data[[symbol_col[1]]])
  symbols <- str_split_fixed(symbols, " /// | // |;|,", 2)[, 1]
  symbols <- str_trim(symbols)
  symbols[symbols == "" | is.na(symbols)] <- NA_character_

  as.data.frame(expr) |>
    mutate(symbol = symbols) |>
    filter(!is.na(symbol)) |>
    group_by(symbol) |>
    summarise(across(where(is.numeric), mean), .groups = "drop") |>
    column_to_rownames("symbol") |>
    as.matrix()
}

fit_external_gene <- function(symbol, gene_expr, clinical, tcga_direction, manual_tier) {
  if (!symbol %in% rownames(gene_expr)) {
    return(tibble(
      symbol = symbol,
      manual_tier = manual_tier,
      external_present = FALSE,
      external_n = nrow(clinical),
      external_events = sum(clinical$os_event == 1, na.rm = TRUE),
      external_log_hr = NA_real_,
      external_hr = NA_real_,
      external_hr_ci_low = NA_real_,
      external_hr_ci_high = NA_real_,
      external_p_value = NA_real_,
      external_adjusted_log_hr = NA_real_,
      external_adjusted_p_value = NA_real_,
      external_same_direction = FALSE
    ))
  }

  dat <- clinical |>
    mutate(
      expr = as.numeric(scale(gene_expr[symbol, sample])),
      grade_high = grade >= 3,
      t_high = t_stage >= 3
    ) |>
    filter(!is.na(os_time), os_time > 0, !is.na(os_event), !is.na(expr))

  fit <- coxph(Surv(os_time, os_event) ~ expr, data = dat)
  term <- tidy(fit, conf.int = TRUE) |> filter(term == "expr")

  adjusted_fit <- tryCatch(
    coxph(Surv(os_time, os_event) ~ expr + age + grade_high + t_high, data = dat),
    error = function(e) NULL
  )
  adjusted_term <- if (is.null(adjusted_fit)) {
    NULL
  } else {
    tidy(adjusted_fit, conf.int = TRUE) |> filter(term == "expr")
  }

  tibble(
    symbol = symbol,
    manual_tier = manual_tier,
    external_present = TRUE,
    external_n = fit$n,
    external_events = fit$nevent,
    external_log_hr = term$estimate,
    external_hr = exp(term$estimate),
    external_hr_ci_low = exp(term$conf.low),
    external_hr_ci_high = exp(term$conf.high),
    external_p_value = term$p.value,
    external_adjusted_log_hr = if (is.null(adjusted_term)) NA_real_ else adjusted_term$estimate,
    external_adjusted_p_value = if (is.null(adjusted_term)) NA_real_ else adjusted_term$p.value,
    external_same_direction = sign(term$estimate) == sign(tcga_direction)
  )
}

gse_path <- file.path(DIRS$processed, "gse29609_series_matrix.rds")
gse <- if (file.exists(gse_path)) {
  readRDS(gse_path)
} else {
  GEOquery::getGEO(
    "GSE29609",
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    destdir = file.path(DIRS$raw, "geo")
  )[[1]]
}

write_rds_atomic(gse, gse_path)

expr <- Biobase::exprs(gse)
pdata <- Biobase::pData(gse)
fdata <- Biobase::fData(gse)
gene_expr <- collapse_gse29609_to_gene(expr, fdata)

clinical <- tibble(
  sample = rownames(pdata),
  os_time = suppressWarnings(as.numeric(pdata[["survival time:ch1"]])),
  os_event = suppressWarnings(as.numeric(pdata[["death (1=yes, 0=no):ch1"]])),
  cancer_death_event = suppressWarnings(as.numeric(pdata[["death from cancer (1=yes, 0=no):ch1"]])),
  age = suppressWarnings(as.numeric(pdata[["age at diagnosis (y):ch1"]])),
  grade = suppressWarnings(as.numeric(pdata[["fuhrman grade:ch1"]])),
  t_stage = suppressWarnings(as.numeric(pdata[["t (tnm stage):ch1"]])),
  n_stage = suppressWarnings(as.numeric(pdata[["n (tnm stage):ch1"]])),
  m_stage = suppressWarnings(as.numeric(pdata[["m (tnm stage):ch1"]]))
)

write_csv_atomic(
  clinical |>
    count(os_event, name = "n_samples") |>
    mutate(accession = "GSE29609", endpoint = "overall survival"),
  file.path(DIRS$tables, "gse29609_sample_summary.csv")
)

candidate_priority <- read_csv(file.path(DIRS$tables, "manuscript_candidate_prioritization.csv"), show_col_types = FALSE)

external_results <- bind_rows(lapply(seq_len(nrow(candidate_priority)), function(i) {
  fit_external_gene(
    symbol = candidate_priority$symbol[i],
    gene_expr = gene_expr,
    clinical = clinical,
    tcga_direction = candidate_priority$main_log_hr[i],
    manual_tier = candidate_priority$manual_tier[i]
  )
})) |>
  mutate(
    external_fdr = p.adjust(external_p_value, method = "BH"),
    external_directional_nominal_support = external_present &
      external_same_direction &
      !is.na(external_p_value) &
      external_p_value < 0.05,
    external_interpretation = case_when(
      !external_present ~ "not present on GSE29609 platform annotation",
      external_directional_nominal_support ~ "same-direction nominal external support",
      external_same_direction ~ "same direction but not nominally significant",
      !external_same_direction & !is.na(external_p_value) & external_p_value < 0.05 ~ "opposite-direction nominal external signal",
      TRUE ~ "not externally supported in this small cohort"
    )
  ) |>
  left_join(
    candidate_priority |>
      select(symbol, main_log_hr, main_hr, main_fdr, manuscript_role, survival_direction, tumor_direction),
    by = "symbol"
  ) |>
  arrange(
    factor(manual_tier, levels = c("lead", "supporting", "supporting risk", "interpret cautiously", "composition flag", "do not highlight")),
    external_p_value
  )

write_csv_atomic(external_results, file.path(DIRS$tables, "external_survival_gse29609.csv"))

external_summary <- tibble(
  metric = c(
    "gse29609_samples",
    "gse29609_events",
    "tested_candidates",
    "platform_present_candidates",
    "same_direction_candidates",
    "same_direction_nominal_candidates",
    "opposite_direction_nominal_candidates"
  ),
  value = c(
    nrow(clinical),
    sum(clinical$os_event == 1, na.rm = TRUE),
    nrow(external_results),
    sum(external_results$external_present, na.rm = TRUE),
    sum(external_results$external_same_direction, na.rm = TRUE),
    sum(external_results$external_directional_nominal_support, na.rm = TRUE),
    sum(!external_results$external_same_direction &
          !is.na(external_results$external_p_value) &
          external_results$external_p_value < 0.05, na.rm = TRUE)
  )
)

write_csv_atomic(external_summary, file.path(DIRS$tables, "external_survival_gse29609_summary.csv"))

message("GSE29609 external survival validation complete.")
