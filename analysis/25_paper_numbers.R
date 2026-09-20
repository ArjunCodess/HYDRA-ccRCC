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
null_summary <- tab("nested_cv_clinical_null_summary.csv")
funnel_test <- tab("funnel_external_paired_bootstrap.csv")
shape <- tab("candidate_survival_shape_sensitivity.csv")
overlap <- tab("null_overlap_check.csv")
checkmate <- tab("checkmate025_study_summary.csv")
cv_conditional <- tab("candidate_cv_clinical_increment.csv")

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
  macro("NestedDelta", sprintf("%.3f", nested$mean_delta_c)),
  macro("NestedCiLow", sprintf("%.3f", nested$patient_bootstrap_ci_low)),
  macro("NestedCiHigh", sprintf("%.3f", nested$patient_bootstrap_ci_high)),
  macro("NestedBrierDelta", sprintf("%.3f", nested$mean_delta_brier3)),
  macro("NestedNoGeneFolds", number(nested$no_gene_folds)),
  macro("NestedOutlierFallbacks", number(nested$outlier_replacement_fallbacks)),
  macro("NestedPassed", ifelse(nested$cv_acceptance, "passed", "failed")),
  macro("NullSelected", number(null_summary$gene_selected)),
  macro("GseFunnelDelta", sprintf("%.3f", funnel_test$difference[funnel_test$cohort == "GSE29609"])),
  macro("GseFunnelCiLow", sprintf("%.3f", funnel_test$ci_low[funnel_test$cohort == "GSE29609"])),
  macro("GseFunnelCiHigh", sprintf("%.3f", funnel_test$ci_high[funnel_test$cohort == "GSE29609"])),
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
  macro("CheckmateCandidates", number(checkmate$value[checkmate$metric == "candidates_mapped"]))
)
writeLines(c("% Generated from committed result tables by analysis/25_paper_numbers.R.", items),
           "paper/results_macros.tex")
message("Wrote ", length(items), " traceable manuscript numbers.")
