source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(readr)
})

required_files <- c(
  file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"),
  file.path(DIRS$tables, "gse40435_sample_summary.csv"),
  file.path(DIRS$tables, "gse53757_sample_summary.csv"),
  file.path(DIRS$tables, "gse29609_sample_summary.csv"),
  FILES$tcga_deg,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"),
  FILES$tcga_survival,
  FILES$tcga_enrichment,
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  file.path(DIRS$tables, "strict_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_ranked_shortlist.csv"),
  file.path(DIRS$tables, "candidate_survival_report.csv"),
  file.path(DIRS$tables, "threshold_sensitivity.csv"),
  file.path(DIRS$tables, "null_overlap_check.csv"),
  file.path(DIRS$tables, "deg_vs_prognostic_comparison.csv"),
  file.path(DIRS$tables, "cell_type_sanity_check.csv"),
  file.path(DIRS$tables, "high_confidence_literature_table.csv"),
  file.path(DIRS$tables, "manuscript_candidate_prioritization.csv"),
  file.path(DIRS$tables, "candidate_clinical_composition_sensitivity.csv"),
  file.path(DIRS$tables, "composition_marker_score_availability.csv"),
  file.path(DIRS$tables, "external_survival_gse29609.csv"),
  file.path(DIRS$tables, "external_survival_gse29609_summary.csv"),
  file.path("results", "high_confidence_gene_dossiers.md"),
  file.path(DIRS$tables, "evidence_funnel.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_pathway_summary.csv"),
  file.path(DIRS$tables, "candidate_summary.csv"),
  file.path(DIRS$figures, "tcga_kirc_pca.png"),
  file.path(DIRS$figures, "tcga_kirc_volcano.png"),
  file.path(DIRS$figures, "tcga_kirc_discordance.png"),
  file.path(DIRS$figures, "tcga_kirc_directional_discordance.png"),
  file.path(DIRS$figures, "candidate_forest_plot.png"),
  file.path(DIRS$figures, "evidence_funnel.png")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required pipeline outputs: ", paste(missing_files, collapse = ", "))
}

if (dir.exists(file.path("results", "logs"))) {
  stop("results/logs exists, but the pipeline must not create a logs directory.")
}

summary <- read_csv(file.path(DIRS$tables, "candidate_summary.csv"), show_col_types = FALSE)
required_metrics <- c(
  "reproducible_deg",
  "main_stage_grade_complete_prognostic",
  "ph_pass",
  "sensitivity_pass",
  "strict_candidate",
  "high_confidence_candidate"
)

missing_metrics <- setdiff(required_metrics, summary$metric)
if (length(missing_metrics) > 0) {
  stop("candidate_summary.csv is missing metrics: ", paste(missing_metrics, collapse = ", "))
}

values <- setNames(summary$value, summary$metric)
if (values[["high_confidence_candidate"]] > values[["strict_candidate"]]) {
  stop("High-confidence candidate count exceeds strict candidate count.")
}
if (values[["strict_candidate"]] > values[["sensitivity_pass"]]) {
  stop("Strict candidate count exceeds sensitivity-pass count.")
}

ranked <- read_csv(file.path(DIRS$tables, "high_confidence_ranked_shortlist.csv"), show_col_types = FALSE)
if (nrow(ranked) != values[["high_confidence_candidate"]]) {
  stop("Ranked shortlist row count does not match high-confidence candidate count.")
}

survival_report <- read_csv(file.path(DIRS$tables, "candidate_survival_report.csv"), show_col_types = FALSE)
required_survival_cols <- c("main_hr", "main_hr_ci_low", "main_hr_ci_high", "main_ph_p_value", "ph_status")
missing_survival_cols <- setdiff(required_survival_cols, names(survival_report))
if (length(missing_survival_cols) > 0) {
  stop("candidate_survival_report.csv is missing columns: ", paste(missing_survival_cols, collapse = ", "))
}

composition <- read_csv(file.path(DIRS$tables, "candidate_clinical_composition_sensitivity.csv"), show_col_types = FALSE)
required_composition_cols <- c(
  "symbol",
  "gene_lrt_fdr_vs_clinical",
  "composition_adjusted_fdr",
  "same_direction_after_composition"
)
missing_composition_cols <- setdiff(required_composition_cols, names(composition))
if (length(missing_composition_cols) > 0) {
  stop("candidate_clinical_composition_sensitivity.csv is missing columns: ", paste(missing_composition_cols, collapse = ", "))
}

external <- read_csv(file.path(DIRS$tables, "external_survival_gse29609.csv"), show_col_types = FALSE)
required_external_cols <- c(
  "symbol",
  "external_present",
  "external_log_hr",
  "external_p_value",
  "external_same_direction",
  "external_interpretation"
)
missing_external_cols <- setdiff(required_external_cols, names(external))
if (length(missing_external_cols) > 0) {
  stop("external_survival_gse29609.csv is missing columns: ", paste(missing_external_cols, collapse = ", "))
}

external_summary <- read_csv(file.path(DIRS$tables, "external_survival_gse29609_summary.csv"), show_col_types = FALSE)
if (!all(c("gse29609_samples", "gse29609_events") %in% external_summary$metric)) {
  stop("external_survival_gse29609_summary.csv is missing required metrics.")
}

message("Output validation complete.")
