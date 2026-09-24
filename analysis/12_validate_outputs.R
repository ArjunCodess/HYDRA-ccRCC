source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(readr)
})

required_files <- c(
  file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"),
  file.path(DIRS$tables, "tcga_kirc_sample_selection_audit.csv"),
  file.path(DIRS$tables, "gse40435_sample_summary.csv"),
  file.path(DIRS$tables, "gse53757_sample_summary.csv"),
  file.path(DIRS$tables, "gse29609_sample_summary.csv"),
  FILES$tcga_deg,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"),
  FILES$tcga_survival,
  FILES$tcga_apeglm_survival,
  FILES$tcga_apeglm_survival_summary,
  FILES$tcga_enrichment,
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  file.path(DIRS$tables, "candidate_paired_de_sensitivity.csv"),
  file.path(DIRS$tables, "prior_candidate_delta.csv"),
  file.path(DIRS$tables, "strict_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_evidence.csv"),
  file.path(DIRS$tables, "candidate_survival_report.csv"),
  file.path(DIRS$tables, "threshold_sensitivity.csv"),
  file.path(DIRS$tables, "null_overlap_check.csv"),
  file.path(DIRS$tables, "deg_vs_prognostic_comparison.csv"),
  file.path(DIRS$tables, "cell_type_sanity_check.csv"),
  file.path(DIRS$tables, "high_confidence_literature_table.csv"),
  file.path(DIRS$tables, "candidate_interpretation_context.csv"),
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
  file.path(DIRS$tables, "nested_cv_predictions.csv"),
  file.path(DIRS$tables, "nested_cv_folds.csv"),
  file.path(DIRS$tables, "nested_cv_repeat_metrics.csv"),
  file.path(DIRS$tables, "nested_cv_patient_bootstrap.csv"),
  file.path(DIRS$tables, "nested_cv_summary.csv"),
  file.path(DIRS$tables, "nested_cv_clinical_null.csv"),
  file.path(DIRS$tables, "nested_cv_clinical_null_summary.csv"),
  file.path(DIRS$tables, "acceptance_criteria.csv"),
  file.path(DIRS$tables, "hpa_candidate_top_cell_types.csv"),
  file.path(DIRS$tables, "hpa_candidate_cell_source_summary.csv"),
  file.path(DIRS$tables, "candidate_direct_tumor_purity_sensitivity.csv"),
  file.path(DIRS$tables, "candidate_survival_shape_sensitivity.csv"),
  file.path(DIRS$tables, "funnel_matched_lists.csv"),
  file.path(DIRS$tables, "funnel_external_gene_results.csv"),
  file.path(DIRS$tables, "funnel_external_summary.csv"),
  file.path(DIRS$tables, "funnel_external_paired_bootstrap.csv"),
  file.path(DIRS$tables, "tumor_purity_coverage.csv"),
  file.path(DIRS$tables, "tracerx_multiregion_source_files.csv"),
  file.path(DIRS$tables, "tracerx_candidate_patient_region_discordance.csv"),
  file.path(DIRS$tables, "tracerx_candidate_multiregion_summary.csv"),
  file.path(DIRS$tables, "tracerx_one_region_cox_repeats.csv"),
  file.path(DIRS$tables, "tracerx_one_region_cox_summary.csv"),
  file.path(DIRS$tables, "tracerx_fixed_subset_patients.csv"),
  file.path(DIRS$tables, "tracerx_multiregion_study_summary.csv"),
  file.path(DIRS$tables, "checkmate025_source_file.csv"),
  file.path(DIRS$tables, "checkmate025_candidate_treatment_interactions.csv"),
  file.path(DIRS$tables, "checkmate025_overall_treatment_effects.csv"),
  file.path(DIRS$tables, "checkmate025_study_summary.csv"),
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
  file.path(DIRS$figures, "evidence_funnel.png"),
  file.path(DIRS$figures, "nested_cv_increment.png"),
  file.path(DIRS$figures, "funnel_external_comparison.png"),
  file.path(DIRS$figures, "nested_selection_benchmark.png"),
  file.path(DIRS$tables, "nested_benchmark_summary.csv"),
  file.path(DIRS$tables, "nested_benchmark_folds.csv"),
  file.path(DIRS$tables, "nested_benchmark_repeat_metrics.csv"),
  file.path(DIRS$tables, "nested_benchmark_predictions.csv"),
  file.path(DIRS$tables, "nested_benchmark_patient_bootstrap.csv"),
  file.path("paper", "results_macros.tex"),
  file.path(DIRS$tables, "gse53757_pair_audit.csv"),
  file.path(DIRS$tables, "gse53757_unpaired_summary.csv"),
  file.path(DIRS$tables, "gse53757_unpaired_sensitivity.csv"),
  file.path(DIRS$tables, "tcga_aliquot_selection_audit.csv"),
  file.path(DIRS$tables, "tcga_aliquot_arithmetic.csv"),
  file.path(DIRS$tables, "aliquot_sensitivity_summary.csv"),
  file.path(DIRS$tables, "published_signature_summary.csv"),
  file.path(DIRS$tables, "clearcode34_symbol_map.csv"),
  file.path(DIRS$tables, "survival_permutation_null_summary.csv"),
  file.path(DIRS$tables, "nested_cv_fold_events.csv"),
  file.path(DIRS$tables, "cohort_dictionary.csv"),
  file.path(DIRS$tables, "clinical_factor_counts.csv"),
  file.path(DIRS$figures, "external_log_hr_forest.png"),
  file.path(DIRS$tables, "gene_mapping_rules.csv"),
  file.path(DIRS$tables, "nested_gene_selection_frequency.csv"),
  file.path(DIRS$tables, "composition_adjustment_display.csv"),
  file.path(DIRS$tables, "tracerx_direction_dispersion.csv"),
  file.path(DIRS$figures, "nested_gene_selection_frequency.png"),
  file.path(DIRS$figures, "composition_marker_fdr.png"),
  file.path(DIRS$figures, "tracerx_sampling_distributions.png"),
  file.path(DIRS$tables, "candidate_evidence_matrix.csv"),
  file.path(DIRS$tables, "ridge_specification.csv"),
  file.path(DIRS$tables, "ridge_fold_audit.csv"),
  file.path(DIRS$tables, "ridge_multigene_summary.csv"),
  file.path(DIRS$tables, "numerical_consistency.csv"),
  file.path(DIRS$figures, "master_funnel_benchmark.png"),
  file.path(DIRS$figures, "gse53757_pairing_sensitivity.png")
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

tcga_deg <- read_csv(FILES$tcga_deg, show_col_types = FALSE)
tcga_samples <- read_csv(file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"), show_col_types = FALSE)
selection_audit <- read_csv(file.path(DIRS$tables, "tcga_kirc_sample_selection_audit.csv"),
                            show_col_types = FALSE)
if (nrow(selection_audit) != 2L ||
    !setequal(selection_audit$shortLetterCode, c("TP", "NT")) ||
    any(selection_audit$selected_samples != selection_audit$unique_patients) ||
    any(selection_audit$removed_replicate_samples !=
        selection_audit$raw_samples - selection_audit$selected_samples) ||
    selection_audit$selected_samples[selection_audit$shortLetterCode == "TP"] != 533L ||
    selection_audit$raw_samples[selection_audit$shortLetterCode == "TP"] != 541L) {
  stop("TCGA sample selection audit is inconsistent with the cached cohort.")
}
if (tcga_samples$n_samples[tcga_samples$sample_type == "Primary Tumor"] != 533L ||
    tcga_samples$n_samples[tcga_samples$sample_type == "Solid Tissue Normal"] != 72L) {
  stop("TCGA selected sample counts differ from the audited cached cohort.")
}
paired_summary <- read_csv(file.path(DIRS$tables, "tcga_kirc_paired_deg_summary.csv"), show_col_types = FALSE)
paired_values <- setNames(paired_summary$value, paired_summary$metric)
if (paired_values[["paired_samples"]] != 2L * paired_values[["paired_patients"]] ||
    paired_values[["paired_patients"]] > 72L) {
  stop("Paired TCGA differential expression has invalid patient coverage.")
}
survival <- read_csv(FILES$tcga_survival, show_col_types = FALSE)
for (model in unique(survival$model_type)) {
  rows <- survival$model_type == model
  if (!isTRUE(all.equal(survival$fdr[rows], p.adjust(survival$p_value[rows], "BH"),
                        tolerance = 1e-12, check.attributes = FALSE))) {
    stop("TCGA survival FDR differs from BH correction for model: ", model)
  }
}
required_apeglm_columns <- c("log2FoldChange_apeglm", "lfcSE_apeglm")
if (!all(required_apeglm_columns %in% names(tcga_deg))) {
  stop("TCGA differential-expression output is missing apeglm MAP estimates.")
}
if (any(!is.finite(tcga_deg$log2FoldChange_apeglm))) {
  stop("TCGA apeglm MAP estimates contain non-finite values.")
}

apeglm_sensitivity <- read_csv(FILES$tcga_apeglm_survival, show_col_types = FALSE)
required_apeglm_sensitivity_columns <- c(
  "gene_id",
  "tcga_mle_log2fc",
  "tcga_apeglm_log2fc",
  "primary_tcga_deg_gate",
  "primary_reproducible_deg",
  "log_hr",
  "p_value",
  "global_fdr",
  "model_status",
  "global_fdr_lt_0_05"
)
missing_apeglm_sensitivity_columns <- setdiff(
  required_apeglm_sensitivity_columns,
  names(apeglm_sensitivity)
)
if (length(missing_apeglm_sensitivity_columns) > 0) {
  stop(
    "All-gene apeglm sensitivity output is missing columns: ",
    paste(missing_apeglm_sensitivity_columns, collapse = ", ")
  )
}
if (anyDuplicated(apeglm_sensitivity$gene_id)) {
  stop("All-gene apeglm sensitivity output contains duplicate gene IDs.")
}
modeled_apeglm <- apeglm_sensitivity$model_status == "ok"
if (any(!is.finite(apeglm_sensitivity$p_value[modeled_apeglm])) ||
    any(!is.finite(apeglm_sensitivity$global_fdr[modeled_apeglm]))) {
  stop("All-gene apeglm sensitivity contains non-finite modeled p-values or FDR values.")
}
expected_global_fdr <- p.adjust(
  apeglm_sensitivity$p_value[modeled_apeglm],
  method = "BH"
)
if (!isTRUE(all.equal(
  apeglm_sensitivity$global_fdr[modeled_apeglm],
  expected_global_fdr,
  tolerance = 1e-12,
  check.attributes = FALSE
))) {
  stop("All-gene sensitivity FDR is not a single BH correction across all modeled genes.")
}

apeglm_summary <- read_csv(FILES$tcga_apeglm_survival_summary, show_col_types = FALSE)
required_apeglm_metrics <- c(
  "qc_filtered_genes",
  "successfully_modeled_genes",
  "global_survival_fdr_lt_0_05"
)
if (!all(required_apeglm_metrics %in% apeglm_summary$metric)) {
  stop("All-gene apeglm sensitivity summary is missing required metrics.")
}
apeglm_summary_values <- setNames(apeglm_summary$value, apeglm_summary$metric)
if (apeglm_summary_values[["qc_filtered_genes"]] != nrow(apeglm_sensitivity) ||
    apeglm_summary_values[["successfully_modeled_genes"]] != sum(modeled_apeglm)) {
  stop("All-gene apeglm sensitivity summary does not match its detail table.")
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
  if (nrow(geo_summary) != 2L || any(geo_summary$n_samples != geo_summary$n_patients)) {
    stop(accession, " parsed tumor-normal pairs do not cover exactly one sample per condition.")
  }
}

candidate_evidence <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_evidence.csv"), show_col_types = FALSE)
if (nrow(candidate_evidence) != values[["high_confidence_candidate"]] ||
    anyDuplicated(candidate_evidence$symbol) ||
    !identical(candidate_evidence$symbol, sort(candidate_evidence$symbol)) ||
    any(c("manual_tier", "final_rank_score", "rank") %in% names(candidate_evidence))) {
  stop("Candidate evidence does not cover each high-confidence gene exactly once.")
}
paired_candidates <- read_csv(file.path(DIRS$tables, "candidate_paired_de_sensitivity.csv"), show_col_types = FALSE)
if (nrow(paired_candidates) != values[["high_confidence_candidate"]] ||
    anyDuplicated(paired_candidates$tcga_gene_id) ||
    any(is.na(paired_candidates$same_direction))) {
  stop("Paired DE sensitivity lacks complete candidate coverage.")
}

survival_report <- read_csv(file.path(DIRS$tables, "candidate_survival_report.csv"), show_col_types = FALSE)
required_survival_cols <- c(
  "main_hr",
  "main_hr_ci_low",
  "main_hr_ci_high",
  "main_ph_p_value",
  "ph_diagnostic_p_ge_0_05"
)
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
  "matched_baseline_log_hr",
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
if (any(!is.finite(purity$matched_baseline_log_hr)) ||
    any(purity$same_direction_after_purity !=
        (sign(purity$gene_log_hr) == sign(purity$matched_baseline_log_hr)))) {
  stop("Direct tumor-purity comparisons are inconsistent with matched baseline fits.")
}

nested_folds <- read_csv(file.path(DIRS$tables, "nested_cv_folds.csv"), show_col_types = FALSE)
nested_scores <- read_csv(file.path(DIRS$tables, "nested_cv_repeat_metrics.csv"), show_col_types = FALSE)
nested_summary <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"), show_col_types = FALSE)
nested_predictions <- read_csv(file.path(DIRS$tables, "nested_cv_predictions.csv"), show_col_types = FALSE)
nested_null <- read_csv(file.path(DIRS$tables, "nested_cv_clinical_null.csv"), show_col_types = FALSE)
nested_boot <- read_csv(file.path(DIRS$tables, "nested_cv_patient_bootstrap.csv"), show_col_types = FALSE)
null_summary <- read_csv(file.path(DIRS$tables, "nested_cv_clinical_null_summary.csv"),
                         show_col_types = FALSE)
if (nrow(nested_folds) != 50L || nrow(nested_scores) != 10L ||
    nrow(nested_summary) != 1L || nrow(nested_null) != 200L ||
    nrow(null_summary) != 1L || null_summary$simulations != nrow(nested_null) ||
    null_summary$gene_selected != sum(!is.na(nested_null$selected_gene)) ||
    any(table(nested_predictions$repeat_id) != nested_summary$n_patients) ||
    anyDuplicated(paste(nested_predictions$repeat_id, nested_predictions$patient_barcode)) ||
    sum(nested_folds$no_gene_selected) != nested_summary$no_gene_folds ||
    sum(nested_folds$outlier_replacement_fallback) !=
      nested_summary$outlier_replacement_fallbacks ||
    any(!is.finite(nested_scores$delta_c)) ||
    any(!is.finite(nested_scores$delta_brier3))) {
  stop("Selection-aware CV or clinical-only null has incomplete coverage.")
}
if (!isTRUE(all.equal(nested_summary$mean_delta_c, mean(nested_scores$delta_c),
                      tolerance = 1e-12, check.attributes = FALSE)) ||
    !isTRUE(all.equal(nested_summary$mean_delta_brier3, mean(nested_scores$delta_brier3),
                      tolerance = 1e-12, check.attributes = FALSE)) ||
    nrow(nested_boot) != 1000L || anyDuplicated(nested_boot$draw) ||
    !isTRUE(all.equal(nested_summary$patient_bootstrap_ci_low,
                      as.numeric(quantile(nested_boot$delta_c, 0.025)),
                      tolerance = 1e-12, check.attributes = FALSE)) ||
    !isTRUE(all.equal(nested_summary$patient_bootstrap_ci_high,
                      as.numeric(quantile(nested_boot$delta_c, 0.975)),
                      tolerance = 1e-12, check.attributes = FALSE)) ||
    !isTRUE(all.equal(null_summary$selection_fraction,
                      mean(!is.na(nested_null$selected_gene)),
                      tolerance = 1e-12, check.attributes = FALSE))) {
  stop("Nested CV summary differs from its repeat-level metrics.")
}
benchmark <- read_csv(file.path(DIRS$tables, "nested_benchmark_summary.csv"), show_col_types = FALSE)
benchmark_folds <- read_csv(file.path(DIRS$tables, "nested_benchmark_folds.csv"), show_col_types = FALSE)
benchmark_repeats <- read_csv(file.path(DIRS$tables, "nested_benchmark_repeat_metrics.csv"),
                              show_col_types = FALSE)
benchmark_predictions <- read_csv(file.path(DIRS$tables, "nested_benchmark_predictions.csv"),
                                  show_col_types = FALSE)
expected_arms <- c("clinical", "hydra", "survival_only", "de_only",
                   "ridge_eligible", "matched_control")
if (!setequal(benchmark$strategy, expected_arms) ||
    any(benchmark$primary_hydra_mismatches != 0L) ||
    nrow(benchmark_folds) != 300L || nrow(benchmark_repeats) != 60L ||
    nrow(benchmark_predictions) != 6L * nested_summary$n_patients * 10L ||
    anyDuplicated(paste(benchmark_predictions$repeat_id, benchmark_predictions$strategy,
                        benchmark_predictions$patient_barcode))) {
  stop("Nested selection benchmark outputs are incomplete or do not match the primary HYDRA genes.")
}

tracerx_discordance <- read_csv(
  file.path(DIRS$tables, "tracerx_candidate_multiregion_summary.csv"),
  show_col_types = FALSE
)
required_tracerx_discordance <- c(
  "symbol",
  "patients_total",
  "multiregion_patients",
  "discordant_multiregion_patients",
  "discordant_multiregion_percent"
)
if (!all(required_tracerx_discordance %in% names(tracerx_discordance))) {
  stop("TRACERx discordance summary is missing required columns.")
}
if (nrow(tracerx_discordance) != values[["high_confidence_candidate"]] ||
    any(tracerx_discordance$discordant_multiregion_percent < 0 |
        tracerx_discordance$discordant_multiregion_percent > 100)) {
  stop("TRACERx discordance summary has invalid candidate coverage or percentages.")
}

tracerx_study <- read_csv(
  file.path(DIRS$tables, "tracerx_multiregion_study_summary.csv"),
  show_col_types = FALSE
)
required_tracerx_metrics <- c(
  "tracerx_primary_regions",
  "tracerx_survival_linked_patients",
  "tracerx_survival_events",
  "tracerx_multiregion_patients",
  "candidates_mapped",
  "region_resampling_repeats",
  "size_matched_patients",
  "size_matched_events"
)
if (!all(required_tracerx_metrics %in% tracerx_study$metric)) {
  stop("TRACERx study summary is missing required metrics.")
}
tracerx_values <- setNames(tracerx_study$value, tracerx_study$metric)
if (tracerx_values[["candidates_mapped"]] != values[["high_confidence_candidate"]] ||
    tracerx_values[["region_resampling_repeats"]] != RESAMPLING$tracerx_region_repeats ||
    tracerx_values[["size_matched_patients"]] != RESAMPLING$tracerx_small_cohort_size) {
  stop("TRACERx study summary does not match the configured analysis.")
}

tracerx_repeats <- read_csv(
  file.path(DIRS$tables, "tracerx_one_region_cox_repeats.csv"),
  show_col_types = FALSE
)
expected_tracerx_rows <- values[["high_confidence_candidate"]] *
  RESAMPLING$tracerx_region_repeats * 3
if (nrow(tracerx_repeats) != expected_tracerx_rows ||
    any(!is.finite(tracerx_repeats$log_hr)) ||
    any(!tracerx_repeats$scenario %in% c("full_cohort", "size_matched_39", "fixed_subset_regions"))) {
  stop("TRACERx one-region repeat table is incomplete or invalid.")
}
if (any(tracerx_repeats$n[tracerx_repeats$scenario == "full_cohort"] !=
        tracerx_values[["tracerx_survival_linked_patients"]]) ||
    any(tracerx_repeats$events[tracerx_repeats$scenario == "full_cohort"] !=
        tracerx_values[["tracerx_survival_events"]]) ||
    any(tracerx_repeats$n[tracerx_repeats$scenario == "size_matched_39"] !=
        tracerx_values[["size_matched_patients"]]) ||
    any(tracerx_repeats$events[tracerx_repeats$scenario == "size_matched_39"] !=
        tracerx_values[["size_matched_events"]]) ||
    any(tracerx_repeats$n[tracerx_repeats$scenario == "fixed_subset_regions"] !=
        tracerx_values[["size_matched_patients"]]) ||
    any(tracerx_repeats$events[tracerx_repeats$scenario == "fixed_subset_regions"] !=
        tracerx_values[["size_matched_events"]])) {
  stop("TRACERx resampling scenarios have inconsistent patient or event counts.")
}

tracerx_resampling <- read_csv(
  file.path(DIRS$tables, "tracerx_one_region_cox_summary.csv"),
  show_col_types = FALSE
)
if (nrow(tracerx_resampling) != values[["high_confidence_candidate"]] * 3 ||
    any(tracerx_resampling$successful_repeats != RESAMPLING$tracerx_region_repeats)) {
  stop("TRACERx one-region summary has incomplete candidate-by-scenario coverage.")
}

tracerx_sources <- read_csv(
  file.path(DIRS$tables, "tracerx_multiregion_source_files.csv"),
  show_col_types = FALSE
)
if (nrow(tracerx_sources) != 3 ||
    any(is.na(tracerx_sources$md5) | tracerx_sources$md5 == "") ||
    any(tracerx_sources$commit != TRACERX_DATA_COMMIT)) {
  stop("TRACERx source inventory is incomplete or not pinned to the configured commit.")
}

checkmate_interactions <- read_csv(
  file.path(DIRS$tables, "checkmate025_candidate_treatment_interactions.csv"),
  show_col_types = FALSE
)
required_checkmate_columns <- c(
  "endpoint", "model", "symbol", "n", "events",
  "interaction_log_hr", "interaction_hr", "interaction_hr_ci_low",
  "interaction_hr_ci_high", "interaction_p_value", "interaction_fdr",
  "interaction_ph_p_value", "nivolumab_log_hr", "everolimus_log_hr"
)
if (!all(required_checkmate_columns %in% names(checkmate_interactions))) {
  stop("CheckMate 025 interaction output is missing required columns.")
}
if (nrow(checkmate_interactions) != values[["high_confidence_candidate"]] * 4 ||
    any(!is.finite(checkmate_interactions$interaction_log_hr)) ||
    any(checkmate_interactions$interaction_p_value < 0 |
        checkmate_interactions$interaction_p_value > 1) ||
    any(checkmate_interactions$interaction_fdr < 0 |
        checkmate_interactions$interaction_fdr > 1) ||
    !setequal(checkmate_interactions$endpoint, c("OS", "PFS")) ||
    !setequal(
      checkmate_interactions$model,
      c("randomized_unadjusted", "age_sex_mskcc_adjusted")
    )) {
  stop("CheckMate 025 interaction output has invalid coverage or estimates.")
}

checkmate_study <- read_csv(
  file.path(DIRS$tables, "checkmate025_study_summary.csv"),
  show_col_types = FALSE
)
required_checkmate_metrics <- c(
  "checkmate025_rna_linked_patients", "nivolumab_patients",
  "everolimus_patients", "os_events", "pfs_events",
  "adjusted_complete_cases", "candidates_mapped"
)
if (!all(required_checkmate_metrics %in% checkmate_study$metric)) {
  stop("CheckMate 025 study summary is missing required metrics.")
}
checkmate_values <- setNames(checkmate_study$value, checkmate_study$metric)
if (checkmate_values[["checkmate025_rna_linked_patients"]] != 250 ||
    checkmate_values[["nivolumab_patients"]] != 120 ||
    checkmate_values[["everolimus_patients"]] != 130 ||
    checkmate_values[["os_events"]] != 191 ||
    checkmate_values[["pfs_events"]] != 222 ||
    checkmate_values[["candidates_mapped"]] != values[["high_confidence_candidate"]]) {
  stop("CheckMate 025 cohort counts do not match the public supplementary data.")
}

checkmate_overall <- read_csv(
  file.path(DIRS$tables, "checkmate025_overall_treatment_effects.csv"),
  show_col_types = FALSE
)
if (nrow(checkmate_overall) != 2 ||
    !setequal(checkmate_overall$endpoint, c("OS", "PFS")) ||
    any(!is.finite(checkmate_overall$nivolumab_vs_everolimus_hr))) {
  stop("CheckMate 025 overall treatment-effect output is incomplete.")
}

checkmate_source <- read_csv(
  file.path(DIRS$tables, "checkmate025_source_file.csv"),
  show_col_types = FALSE
)
if (nrow(checkmate_source) != 1 ||
    tolower(checkmate_source$md5) != CHECKMATE_BRAUN_WORKBOOK_MD5) {
  stop("CheckMate 025 source inventory does not match the pinned workbook.")
}

provenance <- read_csv(file.path(DIRS$tables, "source_provenance.csv"), show_col_types = FALSE)
if (!all(c(
  "TCGA-KIRC",
  "GSE40435",
  "GSE53757",
  "GSE29609",
  "E-MTAB-1980",
  "HPA-v25.1",
  "TRACERx-Renal",
  "CheckMate-025-Braun",
  "Aran-2015-CPE"
) %in%
         provenance$source_id)) {
  stop("Source provenance is incomplete.")
}

manifest <- read_csv(file.path(DIRS$tables, "run_manifest.csv"), show_col_types = FALSE)
if (any(is.na(manifest$md5) | manifest$md5 == "")) {
  stop("Run manifest contains missing checksums.")
}
if (any(!file.exists(manifest$path)) ||
    any(unname(tools::md5sum(manifest$path)) != manifest$md5)) {
  stop("Generated output differs from the run manifest.")
}

funnel_lists <- read_csv(file.path(DIRS$tables, "funnel_matched_lists.csv"), show_col_types = FALSE)
if (dplyr::n_distinct(funnel_lists$rule) != 6L ||
    length(unique(table(funnel_lists$rule))) != 1L ||
    anyDuplicated(paste(funnel_lists$rule, funnel_lists$symbol))) {
  stop("Funnel ablation lists are not equally sized and gene-unique.")
}
funnel_external <- read_csv(file.path(DIRS$tables, "funnel_external_gene_results.csv"), show_col_types = FALSE)
if (nrow(funnel_external) != nrow(funnel_lists) * 2L ||
    any(funnel_external$fdr < 0 | funnel_external$fdr > 1, na.rm = TRUE)) {
  stop("Funnel external evaluation has incomplete coverage or invalid FDR.")
}
acceptance <- read_csv(file.path(DIRS$tables, "acceptance_criteria.csv"), show_col_types = FALSE)
funnel_test <- read_csv(file.path(DIRS$tables, "funnel_external_paired_bootstrap.csv"), show_col_types = FALSE)
expected_cv <- nested_summary$mean_delta_c >= 0.01 &&
  nested_summary$patient_bootstrap_ci_low > 0 && nested_summary$mean_delta_brier3 <= 0
if (!identical(nested_summary$cv_acceptance, expected_cv) ||
    acceptance$status[acceptance$criterion == "selection_aware_cv"] !=
      ifelse(expected_cv, "pass", "fail") ||
    acceptance$status[acceptance$criterion == "external_funnel_vs_survival_only"] !=
      ifelse(all(funnel_test$ci_low > 0), "pass", "fail")) {
  stop("Scientific acceptance rows disagree with the underlying results.")
}

funnel <- read_csv(file.path(DIRS$tables, "evidence_funnel.csv"), show_col_types = FALSE)
if (any(diff(funnel$count) > 0, na.rm = TRUE)) {
  stop("evidence_funnel.csv is not monotonic. Funnel counts must not increase across sequential hardening steps.")
}

if (file.exists(file.path(DIRS$tables, "high_confidence_gene_dossiers.md"))) {
  stop("high_confidence_gene_dossiers.md must be written to results/, not results/tables/.")
}

message("Output validation complete.")
