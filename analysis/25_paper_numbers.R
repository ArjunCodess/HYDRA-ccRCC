source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

read_metric <- function(path, metric) {
  x <- read_csv(path, show_col_types = FALSE)
  value <- x$value[x$metric == metric]
  if (length(value) != 1L || !is.finite(value)) stop("Missing paper metric: ", metric)
  value
}
tab <- function(name) read_csv(file.path(DIRS$tables, name), show_col_types = FALSE)
number <- function(x) formatC(x, format = "d", big.mark = ",")
percent <- function(x) sprintf("%.1f", 100 * x)
macro <- function(name, value) paste0("\\newcommand{\\", name, "}{", value, "}")

sample <- tab("tcga_kirc_sample_summary.csv")
sample_audit <- tab("tcga_kirc_sample_selection_audit.csv")
funnel <- tab("candidate_summary.csv")
funnel_value <- function(metric) funnel$value[funnel$metric == metric]
repro <- tab("reproducibility_summary.csv")
allgene <- tab("tcga_kirc_apeglm_all_gene_survival_summary.csv")
paired <- tab("tcga_kirc_paired_deg_summary.csv")
paired_candidates <- tab("candidate_paired_de_sensitivity.csv")
prior_delta <- tab("prior_candidate_delta.csv")
gse <- tab("external_survival_gse29609_summary.csv")
em <- tab("external_survival_emtab1980_summary.csv")
composition <- tab("candidate_clinical_composition_sensitivity.csv")
bootstrap <- tab("candidate_cox_bootstrap_summary.csv")
tracerx <- tab("tracerx_one_region_cox_summary.csv")
tracerx_cohort <- tab("tracerx_multiregion_study_summary.csv")
nested <- tab("nested_cv_summary.csv")
nested_folds <- tab("nested_cv_folds.csv")
nested_repeats <- tab("nested_cv_repeat_metrics.csv")
null_summary <- tab("nested_cv_clinical_null_summary.csv")
funnel_test <- tab("funnel_external_paired_bootstrap.csv")
funnel_summary <- tab("funnel_external_summary.csv")
funnel_present <- function(rule, cohort) {
  x <- funnel_summary$present[funnel_summary$rule == rule & funnel_summary$cohort == cohort]
  if (length(x) != 1L) stop("Missing matched-list coverage: ", rule, ", ", cohort)
  number(x)
}
cohort_dict <- tab("cohort_dictionary.csv")
tracerx_dispersion <- tab("tracerx_direction_dispersion.csv")
dispersion_value <- function(scenario, column) {
  value <- tracerx_dispersion[[column]][tracerx_dispersion$scenario == scenario]
  if (length(value) != 1L || !is.finite(value)) stop("Missing TRACERx dispersion: ", scenario)
  value
}
shape <- tab("candidate_survival_shape_sensitivity.csv")
overlap <- tab("null_overlap_check.csv")
checkmate <- tab("checkmate025_study_summary.csv")
cv_conditional <- tab("candidate_cv_clinical_increment.csv")
benchmark <- tab("nested_benchmark_summary.csv")
benchmark_folds <- tab("nested_benchmark_folds.csv")
bench <- function(strategy) {
  row <- benchmark[benchmark$strategy == strategy, ]
  if (nrow(row) != 1L) stop("Missing nested-benchmark strategy: ", strategy)
  row
}
dec4 <- function(x) sprintf("%.4f", x)
if (any(benchmark$primary_hydra_mismatches != 0L) || nrow(benchmark) != 6L) {
  stop("Nested benchmark did not reproduce the primary HYDRA genes.")
}
ridge_sizes <- benchmark_folds$n_model_genes[benchmark_folds$strategy == "ridge_eligible"]
hydra_b <- bench("hydra")
survival_b <- bench("survival_only")
de_b <- bench("de_only")
control_b <- bench("matched_control")
ridge_b <- bench("ridge_eligible")
clinical_b <- bench("clinical")
aliquot <- tab("aliquot_sensitivity_summary.csv")
aliquot_n <- function(rule) {
  value <- aliquot$n_high_confidence[aliquot$rule == rule]
  if (length(value) != 1L) stop("Missing aliquot rule: ", rule)
  value
}
unpaired_summary <- tab("gse53757_unpaired_summary.csv")
unpaired_value <- function(metric) {
  value <- unpaired_summary$value[unpaired_summary$metric == metric]
  if (length(value) != 1L) stop("Missing unpaired metric: ", metric)
  value
}
pair_audit <- tab("gse53757_pair_audit_summary.csv")
pair_value <- function(metric) pair_audit$value[pair_audit$metric == metric][1]
signature <- tab("published_signature_summary.csv")
sig_row <- function(cohort, arm) {
  row <- signature[signature$cohort == cohort & signature$arm == arm, ]
  if (nrow(row) != 1L) stop("Missing signature arm: ", cohort, " ", arm)
  row
}
clearcode <- sig_row("TCGA-KIRC", "clinical_plus_clearcode34")
em_subtype <- sig_row("E-MTAB-1980", "clinical_plus_subtype")
em_top <- sig_row("E-MTAB-1980", "clinical_plus_top_gene")
em_panel <- sig_row("E-MTAB-1980", "clinical_plus_signed_panel")
perm_summary <- tab("survival_permutation_null_summary.csv")
fold_event_summary <- tab("nested_cv_fold_event_summary.csv")
exclusions <- tab("nested_cohort_exclusions.csv")
exclusion_n <- function(reason) {
  hit <- exclusions$n[exclusions$reason == reason]
  if (length(hit) != 1L) stop("Missing cohort exclusion: ", reason)
  hit
}
fold_min <- fold_event_summary$value[fold_event_summary$metric == "test_events_min"]
fold_max <- fold_event_summary$value[fold_event_summary$metric == "test_events_max"]
if (length(fold_min) != 1L || fold_min != fold_max) {
  stop("Outer test folds do not share one event count.")
}

items <- c(
  macro("OriginalTumorSamples", number(sample_audit$raw_samples[sample_audit$shortLetterCode == "TP"])),
  macro("TumorPatients", number(sample$n_samples[sample$sample_type == "Primary Tumor"])),
  macro("NormalPatients", number(sample$n_samples[sample$sample_type == "Solid Tissue Normal"])),
  macro("PairedPatients", number(paired$value[paired$metric == "paired_patients"])),
  macro("PairedDeg", number(paired$value[paired$metric == "significant_genes"])),
  macro("PairedCandidatePass", number(sum(paired_candidates$paired_significant, na.rm = TRUE))),
  macro("MappedDeg", number(repro$value[repro$metric == "tcga_significant"])),
  macro("ReproDeg", number(funnel_value("reproducible_deg"))),
  macro("MainSurvival", number(funnel_value("main_stage_grade_complete_prognostic"))),
  macro("SensitivityPass", number(funnel_value("sensitivity_pass"))),
  macro("StrictGenes", number(funnel_value("strict_candidate"))),
  macro("HighGenes", number(funnel_value("high_confidence_candidate"))),
  macro("PriorAdded", number(sum(prior_delta$status == "added"))),
  macro("PriorRemoved", number(sum(prior_delta$status == "removed"))),
  macro("AllGeneCount", number(allgene$value[allgene$metric == "qc_filtered_genes"])),
  macro("AllGeneFdr", number(allgene$value[allgene$metric == "global_survival_fdr_lt_0_05"])),
  macro("NullOverlapP", sprintf("%.3f", overlap$empirical_p_greater_equal)),
  macro("GseMapped", number(gse$value[gse$metric == "platform_present_candidates"])),
  macro("GsePatients", number(gse$value[gse$metric == "gse29609_samples"])),
  macro("GseEvents", number(gse$value[gse$metric == "gse29609_events"])),
  macro("GseSame", number(gse$value[gse$metric == "same_direction_candidates"])),
  macro("GseOppositeNominal", number(gse$value[gse$metric == "opposite_direction_nominal_candidates"])),
  macro("EmMapped", number(em$value[em$metric == "platform_present_candidates"])),
  macro("EmPatients", number(em$value[em$metric == "emtab1980_samples"])),
  macro("EmEvents", number(em$value[em$metric == "emtab1980_events"])),
  macro("EmSame", number(em$value[em$metric == "same_direction_candidates"])),
  macro("EmStrict", number(em$value[em$metric == "strict_external_support_candidates"])),
  macro("CompositionFailures", number(sum(composition$composition_adjusted_fdr >= 0.05, na.rm = TRUE))),
  macro("BootstrapIntervals", number(sum(bootstrap$ci_excludes_zero))),
  macro("ConditionalCvPositive", number(sum(cv_conditional$mean_delta_c_index > 0))),
  macro("NestedDelta", sprintf("%.4f", nested$mean_delta_c)),
  macro("NestedCiLow", sprintf("%.4f", nested$patient_bootstrap_ci_low)),
  macro("NestedCiHigh", sprintf("%.4f", nested$patient_bootstrap_ci_high)),
  macro("NestedBrierDelta", sprintf("%.4f", nested$mean_delta_brier3)),
  macro("NestedClinicalC", sprintf("%.3f", mean(nested_repeats$clinical_c))),
  macro("NestedGeneC", sprintf("%.3f", mean(nested_repeats$gene_c))),
  macro("NestedClinicalSlope", sprintf("%.2f", mean(nested_repeats$clinical_calibration_slope))),
  macro("NestedGeneSlope", sprintf("%.2f", mean(nested_repeats$gene_calibration_slope))),
  macro("NestedDistinctGenes", number(n_distinct(nested_folds$selected_gene, na.rm = TRUE))),
  macro("NestedTopGeneFolds", number(max(table(nested_folds$selected_gene)))),
  macro("NestedNoGeneFolds", number(nested$no_gene_folds)),
  macro("NestedOutlierFallbacks", number(nested$outlier_replacement_fallbacks)),
  macro("NestedPassed", ifelse(nested$cv_acceptance, "passed", "failed")),
  macro("NullSelected", number(null_summary$gene_selected)),
  macro("GseFunnelDelta", sprintf("%.3f", funnel_test$difference[funnel_test$cohort == "GSE29609"])),
  macro("GseFunnelCiLow", sprintf("%.3f", funnel_test$ci_low[funnel_test$cohort == "GSE29609"])),
  macro("GseFunnelCiHigh", sprintf("%.3f", funnel_test$ci_high[funnel_test$cohort == "GSE29609"])),
  macro("GseCompletePresent", funnel_present("complete_rule", "GSE29609")),
  macro("GseSurvivalPresent", funnel_present("survival_only", "GSE29609")),
  macro("EmCompletePresent", funnel_present("complete_rule", "E-MTAB-1980")),
  macro("EmSurvivalPresent", funnel_present("survival_only", "E-MTAB-1980")),
  macro("FunnelPassed", ifelse(all(funnel_test$acceptance), "passed", "failed")),
  macro("ShapeNonlinear", number(sum(shape$nonlinear_fdr < 0.05, na.rm = TRUE))),
  macro("ShapeTimeVarying", number(sum(shape$time_varying_fdr < 0.05, na.rm = TRUE))),
  macro("TracerxFixedMedian", percent(median(tracerx$same_tcga_direction_fraction[tracerx$scenario == "fixed_subset_regions"]))),
  macro("TracerxChangingMedian", percent(median(tracerx$same_tcga_direction_fraction[tracerx$scenario == "size_matched_39"]))),
  macro("TracerxRegions", number(tracerx_cohort$value[tracerx_cohort$metric == "tracerx_primary_regions"])),
  macro("TracerxPatients", number(tracerx_cohort$value[tracerx_cohort$metric == "tracerx_survival_linked_patients"])),
  macro("TracerxEvents", number(tracerx_cohort$value[tracerx_cohort$metric == "tracerx_survival_events"])),
  macro("TracerxSubsetPatients", number(tracerx_cohort$value[tracerx_cohort$metric == "size_matched_patients"])),
  macro("TracerxSubsetEvents", number(tracerx_cohort$value[tracerx_cohort$metric == "size_matched_events"])),
  macro("CheckmatePatients", number(checkmate$value[checkmate$metric == "checkmate025_rna_linked_patients"])),
  macro("CheckmateNivolumab", number(checkmate$value[checkmate$metric == "nivolumab_patients"])),
  macro("CheckmateEverolimus", number(checkmate$value[checkmate$metric == "everolimus_patients"])),
  macro("CheckmateCandidates", number(checkmate$value[checkmate$metric == "candidates_mapped"])),
  macro("BenchHydraDelta", dec4(hydra_b$mean_delta_c)),
  macro("BenchHydraCiLow", dec4(hydra_b$patient_bootstrap_ci_low)),
  macro("BenchHydraCiHigh", dec4(hydra_b$patient_bootstrap_ci_high)),
  macro("BenchHydraSlope", sprintf("%.2f", hydra_b$mean_calibration_slope)),
  macro("BenchHydraGenes", number(hydra_b$distinct_genes)),
  macro("BenchSurvivalDelta", dec4(survival_b$mean_delta_c)),
  macro("BenchSurvivalCiLow", dec4(survival_b$patient_bootstrap_ci_low)),
  macro("BenchSurvivalCiHigh", dec4(survival_b$patient_bootstrap_ci_high)),
  macro("BenchSurvivalGenes", number(survival_b$distinct_genes)),
  macro("BenchDeDelta", dec4(de_b$mean_delta_c)),
  macro("BenchDeCiLow", dec4(de_b$patient_bootstrap_ci_low)),
  macro("BenchDeCiHigh", dec4(de_b$patient_bootstrap_ci_high)),
  macro("BenchDeGenes", number(de_b$distinct_genes)),
  macro("BenchControlDelta", dec4(control_b$mean_delta_c)),
  macro("BenchControlCiLow", dec4(control_b$patient_bootstrap_ci_low)),
  macro("BenchControlCiHigh", dec4(control_b$patient_bootstrap_ci_high)),
  macro("BenchControlGenes", number(control_b$distinct_genes)),
  macro("BenchRidgeDelta", dec4(ridge_b$mean_delta_c)),
  macro("BenchRidgeCiLow", dec4(ridge_b$patient_bootstrap_ci_low)),
  macro("BenchRidgeCiHigh", dec4(ridge_b$patient_bootstrap_ci_high)),
  macro("BenchRidgeBrier", dec4(ridge_b$mean_delta_brier3)),
  macro("BenchRidgeBrierCiLow", dec4(ridge_b$patient_bootstrap_brier_ci_low)),
  macro("BenchRidgeBrierCiHigh", dec4(ridge_b$patient_bootstrap_brier_ci_high)),
  macro("BenchRidgeSlope", sprintf("%.2f", ridge_b$mean_calibration_slope)),
  macro("BenchRidgeGenes", number(round(mean(ridge_sizes)))),
  macro("BenchRidgeGenesMin", number(min(ridge_sizes))),
  macro("BenchRidgeGenesMax", number(max(ridge_sizes))),
  macro("RidgeLambdaMedian", sprintf("%.4g", {
    ridge_audit <- tab("ridge_multigene_summary.csv")
    ridge_audit$value[ridge_audit$metric == "lambda_median"]
  })),
  macro("BenchClinicalC", sprintf("%.3f", clinical_b$mean_c)),
  macro("BenchClinicalSlope", sprintf("%.2f", clinical_b$mean_calibration_slope)),
  macro("AliquotDrop", number(aliquot_n("drop_multi"))),
  macro("AliquotFirst", number(aliquot_n("first_barcode"))),
  macro("AliquotSum", number(aliquot_n("sum_counts"))),
  macro("UnpairedHigh", number(unpaired_value("high_confidence_still_reproducible"))),
  macro("UnpairedJaccard", sprintf("%.3f", unpaired_value("reproducible_jaccard"))),
  macro("PairSameToken", number(pair_value("same_title_token_pairs"))),
  macro("PairSameStage", number(pair_value("same_stage_pairs"))),
  macro("ClearcodeDelta", dec4(clearcode$delta_c)),
  macro("ClearcodeCiLow", dec4(clearcode$delta_c_ci_low)),
  macro("ClearcodeCiHigh", dec4(clearcode$delta_c_ci_high)),
  macro("ClearcodeBrier", dec4(clearcode$delta_brier3)),
  macro("ClearcodeSlope", sprintf("%.2f", clearcode$calibration_slope)),
  macro("ClearcodeGenes", number(clearcode$genes_used)),
  macro("EmSubtypeDelta", dec4(em_subtype$delta_c)),
  macro("EmSubtypeCiLow", dec4(em_subtype$delta_c_ci_low)),
  macro("EmSubtypeCiHigh", dec4(em_subtype$delta_c_ci_high)),
  macro("EmTopDelta", dec4(em_top$delta_c)),
  macro("EmPanelDelta", dec4(em_panel$delta_c)),
  macro("EmPanelGenes", number(em_panel$genes_used)),
  macro("PermSelected", number(perm_summary$value[perm_summary$metric == "gene_selected"])),
  macro("FoldTestEvents", number(fold_event_summary$value[fold_event_summary$metric == "test_events_min"])),
  macro("TracerxFixedIqrLow", percent(dispersion_value("fixed_subset_regions", "iqr_low"))),
  macro("TracerxFixedIqrHigh", percent(dispersion_value("fixed_subset_regions", "iqr_high"))),
  macro("TracerxChangingIqrLow", percent(dispersion_value("size_matched_39", "iqr_low"))),
  macro("TracerxChangingIqrHigh", percent(dispersion_value("size_matched_39", "iqr_high"))),
  macro("CompositionFailedSymbols", paste(sort(composition$symbol[composition$composition_adjusted_fdr >= 0.05]), collapse = ", ")),
  macro("TcgaTumorEvents", number(cohort_dict$events[cohort_dict$dataset == "TCGA-KIRC tumors"])),
  macro("NestedPatients", number(exclusion_n("complete_case"))),
  macro("NestedEvents", number(exclusions$deaths[exclusions$reason == "complete_case"])),
  macro("NestedExcluded", number(sum(exclusions$n[exclusions$reason != "complete_case"]))),
  macro("ExcludeGrade", number(exclusion_n("missing_grade"))),
  macro("ExcludeStage", number(exclusion_n("missing_stage"))),
  macro("ExcludeAge", number(exclusion_n("missing_age"))),
  macro("ExcludeTime", number(exclusion_n("zero_or_nonpositive_survival_time"))),
  macro("ExcludeDeaths", number(sum(exclusions$deaths[exclusions$reason != "complete_case"]))),
  macro("StageOnlyPatients", number(exclusion_n("complete_case") + exclusion_n("missing_grade"))),
  macro("GradeOnlyPatients", number(exclusion_n("complete_case") + exclusion_n("missing_stage")))
)
writeLines(c("% Generated from committed result tables by analysis/25_paper_numbers.R.", items),
           "paper/results_macros.tex")
message("Wrote ", length(items), " traceable manuscript numbers.")
