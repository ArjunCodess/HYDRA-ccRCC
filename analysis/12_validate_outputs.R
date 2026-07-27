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
  file.path(DIRS$tables, "external_survival_emtab1980.csv"),
  file.path(DIRS$tables, "external_survival_emtab1980_summary.csv"),
  file.path(DIRS$tables, "candidate_cox_bootstrap_repeats.csv"),
  file.path(DIRS$tables, "candidate_cox_bootstrap_summary.csv"),
  file.path(DIRS$tables, "candidate_cv_clinical_increment.csv"),
  file.path(DIRS$tables, "candidate_cv_clinical_increment_repeats.csv"),
  file.path(DIRS$tables, "hpa_candidate_top_cell_types.csv"),
  file.path(DIRS$tables, "hpa_candidate_cell_source_summary.csv"),
  file.path(DIRS$tables, "candidate_direct_tumor_purity_sensitivity.csv"),
  file.path(DIRS$tables, "tumor_purity_coverage.csv"),
  file.path(DIRS$tables, "source_provenance.csv"),
  file.path(DIRS$tables, "run_manifest.csv"),
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

for (accession in c("gse40435", "gse53757")) {
  geo_result <- read_csv(
    file.path(DIRS$tables, paste0(accession, "_limma_tumor_vs_normal.csv")),
    show_col_types = FALSE
  )
  required_geo_columns <- c("log2fc_ci_low", "log2fc_ci_high")
  if (!all(required_geo_columns %in% names(geo_result))) {
    stop(accession, " differential-expression output is missing confidence intervals.")
  }

  geo_summary <- read_csv(
    file.path(DIRS$tables, paste0(accession, "_sample_summary.csv")),
    show_col_types = FALSE
  )
  required_sva_columns <- c(
    "n_patients",
    "n_surrogate_variables",
    "full_design_rank",
    "adjusted_design_rank"
  )
  if (!all(required_sva_columns %in% names(geo_summary))) {
    stop(accession, " sample summary is missing SVA design diagnostics.")
  }
  if (any(geo_summary$n_surrogate_variables < 0) ||
      any(geo_summary$adjusted_design_rank < geo_summary$full_design_rank)) {
    stop(accession, " contains invalid SVA design diagnostics.")
  }
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
  "gene_log_hr",
  "gene_p_value",
  "gene_fdr",
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

emtab_external <- read_csv(file.path(DIRS$tables, "external_survival_emtab1980.csv"), show_col_types = FALSE)
required_emtab_external_cols <- c(
  "symbol",
  "external_present",
  "external_log_hr",
  "external_p_value",
  "external_ph_p_value",
  "external_adjusted_log_hr",
  "external_adjusted_p_value",
  "external_adjusted_ph_p_value",
  "external_same_direction",
  "external_strict_support",
  "external_interpretation"
)
missing_emtab_external_cols <- setdiff(required_emtab_external_cols, names(emtab_external))
if (length(missing_emtab_external_cols) > 0) {
  stop(
    "external_survival_emtab1980.csv is missing columns: ",
    paste(missing_emtab_external_cols, collapse = ", ")
  )
}

emtab_summary <- read_csv(
  file.path(DIRS$tables, "external_survival_emtab1980_summary.csv"),
  show_col_types = FALSE
)
if (!all(c("emtab1980_samples", "emtab1980_events") %in% emtab_summary$metric)) {
  stop("external_survival_emtab1980_summary.csv is missing required metrics.")
}
emtab_values <- setNames(emtab_summary$value, emtab_summary$metric)
if (emtab_values[["emtab1980_samples"]] != 101) {
  stop("E-MTAB-1980 sample count must be 101.")
}
if (emtab_values[["emtab1980_events"]] != 23) {
  stop("E-MTAB-1980 event count must be 23.")
}

bootstrap_summary <- read_csv(
  file.path(DIRS$tables, "candidate_cox_bootstrap_summary.csv"),
  show_col_types = FALSE
)
if (nrow(bootstrap_summary) != values[["high_confidence_candidate"]]) {
  stop("Cox bootstrap summary count does not match high-confidence candidate count.")
}
required_bootstrap_columns <- c(
  "bootstrap_se",
  "bootstrap_ci_low",
  "bootstrap_ci_high",
  "bootstrap_bias",
  "direction_agreement"
)
if (!all(required_bootstrap_columns %in% names(bootstrap_summary))) {
  stop("Cox bootstrap summary is missing uncertainty columns.")
}
if (any(
  bootstrap_summary$successful_repeats != RESAMPLING$coefficient_bootstrap_repeats |
    !is.finite(bootstrap_summary$bootstrap_se) |
    bootstrap_summary$bootstrap_ci_low > bootstrap_summary$bootstrap_ci_high
)) {
  stop("Cox bootstrap summary contains incomplete or invalid estimates.")
}

bootstrap_repeats <- read_csv(
  file.path(DIRS$tables, "candidate_cox_bootstrap_repeats.csv"),
  show_col_types = FALSE
)
expected_bootstrap_rows <- values[["high_confidence_candidate"]] *
  RESAMPLING$coefficient_bootstrap_repeats
if (nrow(bootstrap_repeats) != expected_bootstrap_rows) {
  stop("Cox bootstrap repeat table does not have complete candidate-by-repeat coverage.")
}

cv <- read_csv(file.path(DIRS$tables, "candidate_cv_clinical_increment.csv"), show_col_types = FALSE)
if (nrow(cv) != values[["high_confidence_candidate"]]) {
  stop("Clinical-increment candidate count does not match high-confidence candidate count.")
}
if (any(cv$mean_clinical_c_index < 0 | cv$mean_clinical_c_index > 1 |
        cv$mean_clinical_gene_c_index < 0 | cv$mean_clinical_gene_c_index > 1, na.rm = TRUE)) {
  stop("Cross-validated concordance values must be between zero and one.")
}

hpa <- read_csv(file.path(DIRS$tables, "hpa_candidate_cell_source_summary.csv"), show_col_types = FALSE)
if (length(unique(hpa$symbol)) != values[["high_confidence_candidate"]]) {
  stop("HPA cell-source output does not cover every high-confidence candidate.")
}

purity <- read_csv(
  file.path(DIRS$tables, "candidate_direct_tumor_purity_sensitivity.csv"),
  show_col_types = FALSE
)
required_purity_columns <- c(
  "symbol",
  "gene_log_hr",
  "gene_p_value",
  "gene_fdr",
  "gene_ph_p_value",
  "same_direction_after_purity",
  "relative_log_hr_attenuation",
  "n",
  "events"
)
missing_purity_columns <- setdiff(required_purity_columns, names(purity))
if (length(missing_purity_columns) > 0) {
  stop(
    "Direct tumor-purity output is missing columns: ",
    paste(missing_purity_columns, collapse = ", ")
  )
}
if (nrow(purity) != values[["high_confidence_candidate"]]) {
  stop("Direct tumor-purity sensitivity must include every revised candidate.")
}
if (any(!is.finite(purity$gene_log_hr) | !is.finite(purity$gene_p_value))) {
  stop("Direct tumor-purity output contains non-finite gene estimates.")
}

provenance <- read_csv(file.path(DIRS$tables, "source_provenance.csv"), show_col_types = FALSE)
if (!all(c(
  "TCGA-KIRC",
  "GSE40435",
  "GSE53757",
  "GSE29609",
  "E-MTAB-1980",
  "HPA-v25.1",
  "Aran-2015-CPE"
) %in%
         provenance$source_id)) {
  stop("Source provenance is incomplete.")
}

manifest <- read_csv(file.path(DIRS$tables, "run_manifest.csv"), show_col_types = FALSE)
if (any(is.na(manifest$md5) | manifest$md5 == "")) {
  stop("Run manifest contains missing checksums.")
}

funnel <- read_csv(file.path(DIRS$tables, "evidence_funnel.csv"), show_col_types = FALSE)
if (any(diff(funnel$count) > 0, na.rm = TRUE)) {
  stop("evidence_funnel.csv is not monotonic. Funnel counts must not increase across sequential hardening steps.")
}

if (file.exists(file.path(DIRS$tables, "high_confidence_gene_dossiers.md"))) {
  stop("high_confidence_gene_dossiers.md must be written to results/, not results/tables/.")
}

message("Output validation complete.")
