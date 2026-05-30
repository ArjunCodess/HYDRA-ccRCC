source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

repro <- read_csv(file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"), show_col_types = FALSE)
surv <- read_csv(FILES$tcga_survival, show_col_types = FALSE)

model_wide <- surv |>
  select(
    gene_id,
    model_type,
    log_hr,
    hr,
    log_hr_ci_low,
    log_hr_ci_high,
    hr_ci_low,
    hr_ci_high,
    std_error,
    p_value,
    fdr,
    ph_p_value,
    n,
    events,
    covariates,
    warning_count,
    warning_text
  ) |>
  pivot_wider(
    names_from = model_type,
    values_from = c(
      log_hr,
      hr,
      log_hr_ci_low,
      log_hr_ci_high,
      hr_ci_low,
      hr_ci_high,
      std_error,
      p_value,
      fdr,
      ph_p_value,
      n,
      events,
      covariates,
      warning_count,
      warning_text
    ),
    names_glue = "{model_type}_{.value}"
  )

safe_neglog10 <- function(x) {
  out <- -log10(pmax(x, .Machine$double.xmin, na.rm = TRUE))
  out[!is.finite(out)] <- NA_real_
  out
}

load_pathway_classes <- function(symbols) {
  empty <- tibble(symbol = symbols, pathway_class = "Unclassified")
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(empty)

  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) NULL
  )
  if (is.null(hallmark) || !all(c("gene_symbol", "gs_name") %in% names(hallmark))) return(empty)

  pathway_map <- hallmark |>
    transmute(
      symbol = gene_symbol,
      pathway_class = case_when(
        gs_name == "HALLMARK_HYPOXIA" ~ "Hypoxia",
        gs_name == "HALLMARK_ANGIOGENESIS" ~ "Angiogenesis",
        gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ "ECM/EMT",
        gs_name %in% c(
          "HALLMARK_GLYCOLYSIS",
          "HALLMARK_FATTY_ACID_METABOLISM",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_ADIPOGENESIS",
          "HALLMARK_BILE_ACID_METABOLISM",
          "HALLMARK_XENOBIOTIC_METABOLISM"
        ) ~ "Metabolism",
        gs_name %in% c(
          "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING",
          "HALLMARK_COMPLEMENT",
          "HALLMARK_ALLOGRAFT_REJECTION",
          "HALLMARK_IL2_STAT5_SIGNALING",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        ) ~ "Immune",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(pathway_class), symbol %in% symbols) |>
    distinct(symbol, pathway_class)

  pathway_map |>
    group_by(symbol) |>
    summarise(pathway_class = paste(sort(unique(pathway_class)), collapse = ";"), .groups = "drop") |>
    right_join(tibble(symbol = symbols), by = "symbol") |>
    mutate(pathway_class = if_else(is.na(pathway_class), "Unclassified", pathway_class))
}

pathway_classes <- load_pathway_classes(unique(repro$symbol))

candidates <- repro |>
  left_join(model_wide, by = c("tcga_gene_id" = "gene_id")) |>
  left_join(pathway_classes, by = "symbol") |>
  mutate(
    main_log_hr = stage_grade_complete_log_hr,
    main_hr = stage_grade_complete_hr,
    main_hr_ci_low = stage_grade_complete_hr_ci_low,
    main_hr_ci_high = stage_grade_complete_hr_ci_high,
    main_log_hr_ci_low = stage_grade_complete_log_hr_ci_low,
    main_log_hr_ci_high = stage_grade_complete_log_hr_ci_high,
    main_std_error = stage_grade_complete_std_error,
    main_p_value = stage_grade_complete_p_value,
    main_fdr = stage_grade_complete_fdr,
    main_ph_p_value = stage_grade_complete_ph_p_value,
    main_n = stage_grade_complete_n,
    main_events = stage_grade_complete_events,
    main_warning_count = stage_grade_complete_warning_count,
    main_warning_text = stage_grade_complete_warning_text,
    stage_sensitivity_same_direction = !is.na(stage_complete_log_hr) & sign(stage_complete_log_hr) == sign(main_log_hr),
    grade_sensitivity_same_direction = !is.na(grade_complete_log_hr) & sign(grade_complete_log_hr) == sign(main_log_hr),
    stage_sensitivity_nominal = !is.na(stage_complete_p_value) & stage_complete_p_value < 0.05,
    grade_sensitivity_nominal = !is.na(grade_complete_p_value) & grade_complete_p_value < 0.05,
    ph_pass = !is.na(main_ph_p_value) & main_ph_p_value >= THRESHOLDS$ph_min_p,
    meaningful_survival_effect = !is.na(main_log_hr) & abs(main_log_hr) >= THRESHOLDS$min_abs_log_hr,
    geo_effect_support = abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc &
      abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
    prognostic = !is.na(main_fdr) & main_fdr < THRESHOLDS$strict_survival_fdr,
    sensitivity_pass = stage_sensitivity_same_direction & grade_sensitivity_same_direction &
      stage_sensitivity_nominal & grade_sensitivity_nominal,
    strict_candidate = reproducible_deg & prognostic & ph_pass & meaningful_survival_effect &
      geo_effect_support & sensitivity_pass,
    high_confidence_candidate = strict_candidate &
      main_fdr < THRESHOLDS$high_confidence_survival_fdr &
      abs(main_log_hr) >= THRESHOLDS$high_confidence_abs_log_hr,
    candidate = strict_candidate,
    abs_tcga_log2fc = abs(tcga_log2fc),
    abs_log_hr = abs(main_log_hr),
    evidence_score = safe_neglog10(main_fdr) +
      abs(main_log_hr) +
      pmin(abs(tcga_log2fc), 5) / 5 +
      pmin(abs(gse40435_log2fc), 3) / 3 +
      pmin(abs(gse53757_log2fc), 3) / 3
  ) |>
  arrange(desc(strict_candidate), desc(evidence_score), main_fdr)

write_csv_atomic(candidates, file.path(DIRS$tables, "candidate_gene_evidence_table.csv"))

strict_candidates <- candidates |>
  filter(strict_candidate) |>
  select(
    symbol,
    pathway_class,
    evidence_score,
    tcga_gene_id,
    tcga_log2fc,
    gse40435_log2fc,
    gse53757_log2fc,
    main_log_hr,
    main_hr,
    main_hr_ci_low,
    main_hr_ci_high,
    main_fdr,
    main_ph_p_value,
    main_warning_count,
    stage_complete_log_hr,
    stage_complete_p_value,
    grade_complete_log_hr,
    grade_complete_p_value
  ) |>
  arrange(desc(evidence_score))

write_csv_atomic(strict_candidates, file.path(DIRS$tables, "strict_candidate_genes.csv"))

high_confidence_candidates <- candidates |>
  filter(high_confidence_candidate) |>
  select(
    symbol,
    pathway_class,
    evidence_score,
    tcga_gene_id,
    tcga_log2fc,
    gse40435_log2fc,
    gse53757_log2fc,
    main_log_hr,
    main_hr,
    main_hr_ci_low,
    main_hr_ci_high,
    main_fdr,
    main_ph_p_value,
    main_warning_count,
    stage_complete_log_hr,
    stage_complete_p_value,
    grade_complete_log_hr,
    grade_complete_p_value
  ) |>
  arrange(desc(evidence_score))

write_csv_atomic(high_confidence_candidates, file.path(DIRS$tables, "high_confidence_candidate_genes.csv"))

pathway_summary <- candidates |>
  filter(high_confidence_candidate) |>
  tidyr::separate_rows(pathway_class, sep = ";") |>
  count(pathway_class, sort = TRUE, name = "n_high_confidence_candidates")

write_csv_atomic(pathway_summary, file.path(DIRS$tables, "high_confidence_candidate_pathway_summary.csv"))

summary <- tibble(
  metric = c(
    "reproducible_deg",
    "main_stage_grade_complete_prognostic",
    "ph_pass",
    "sensitivity_pass",
    "strict_candidate",
    "high_confidence_candidate"
  ),
  value = c(
    sum(candidates$reproducible_deg, na.rm = TRUE),
    sum(candidates$reproducible_deg & candidates$prognostic, na.rm = TRUE),
    sum(candidates$reproducible_deg & candidates$prognostic & candidates$ph_pass, na.rm = TRUE),
    sum(candidates$reproducible_deg & candidates$prognostic & candidates$ph_pass & candidates$sensitivity_pass, na.rm = TRUE),
    sum(candidates$strict_candidate, na.rm = TRUE),
    sum(candidates$high_confidence_candidate, na.rm = TRUE)
  )
)
write_csv_atomic(summary, file.path(DIRS$tables, "candidate_summary.csv"))

message("Strict candidate genes: ", sum(candidates$strict_candidate, na.rm = TRUE))
message("High-confidence candidate genes: ", sum(candidates$high_confidence_candidate, na.rm = TRUE))
