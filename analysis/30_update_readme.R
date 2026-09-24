source("analysis/00_config.R")
suppressPackageStartupMessages(library(readr))

nested <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"),
                   show_col_types = FALSE)
null <- read_csv(file.path(DIRS$tables, "nested_cv_clinical_null_summary.csv"),
                 show_col_types = FALSE)
folds <- read_csv(file.path(DIRS$tables, "nested_cv_folds.csv"),
                  show_col_types = FALSE)
repeats <- read_csv(file.path(DIRS$tables, "nested_cv_repeat_metrics.csv"),
                    show_col_types = FALSE)
if (nrow(nested) != 1L || nested$repeats != 10L || nested$folds != 5L ||
    nrow(null) != 1L || null$simulations != 200L ||
    nrow(folds) != 50L || nrow(repeats) != 10L) {
  stop("Final nested CV and null summaries are required before updating README.")
}

line <- paste0(
  "The [selection-aware nested CV](results/tables/nested_cv_summary.csv) used ",
  nested$n_patients, " patients across ten repeats of five folds. It selected no gene in ",
  nested$no_gene_folds, " folds; ", nested$outlier_replacement_fallbacks,
  " folds required a DESeq2 no-replacement retry. The top-ranked gene varied ",
  "across ", dplyr::n_distinct(folds$selected_gene, na.rm = TRUE),
  " genes, with the most frequent selected in ", max(table(folds$selected_gene)),
  " folds. The [selection-frequency figure](results/figures/nested_gene_selection_frequency.png) shows those ten genes. The mean selected-gene minus ",
  "clinical concordance was ", sprintf("%+.4f", nested$mean_delta_c),
  " (patient-resampled 95% interval ",
  sprintf("%+.4f", nested$patient_bootstrap_ci_low), " to ",
  sprintf("%+.4f", nested$patient_bootstrap_ci_high),
  "), and the three-year Brier-score difference was ",
  sprintf("%+.4f", nested$mean_delta_brier3), ". The prespecified CV criterion **",
  ifelse(nested$cv_acceptance, "passed", "failed"),
  "**. Mean risk-score calibration slopes were ",
  sprintf("%.2f", mean(repeats$clinical_calibration_slope)), " for clinical-only and ",
  sprintf("%.2f", mean(repeats$gene_calibration_slope)),
  " for the selected-gene models, so discrimination and calibration did not move together. ",
  "The [200 clinical-only null simulations](results/tables/nested_cv_clinical_null_summary.csv) ",
  "selected a gene in ", null$gene_selected,
  " single-fold runs under their simulated clinical-risk model. The patient ",
  "bootstrap keeps the fitted fold models fixed, ",
  "so its interval omits training-set and split uncertainty."
)

path <- "README.md"
lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
start <- which(lines == "<!-- nested-results:start -->")
end <- which(lines == "<!-- nested-results:end -->")
if (length(start) != 1L || length(end) != 1L || end <= start) {
  stop("README nested-result markers are missing or out of order.")
}
writeLines(c(lines[seq_len(start)], line, lines[end:length(lines)]),
           path, useBytes = TRUE)

benchmark <- read_csv(file.path(DIRS$tables, "nested_benchmark_summary.csv"),
                      show_col_types = FALSE)
benchmark_folds <- read_csv(file.path(DIRS$tables, "nested_benchmark_folds.csv"),
                            show_col_types = FALSE)
arm <- function(strategy) {
  row <- benchmark[benchmark$strategy == strategy, ]
  if (nrow(row) != 1L) stop("Missing benchmark strategy: ", strategy)
  row
}
fmt <- function(x) sprintf("%+.4f", x)
interval <- function(row) {
  paste0(fmt(row$patient_bootstrap_ci_low), " to ", fmt(row$patient_bootstrap_ci_high))
}
hydra <- arm("hydra")
survival_arm <- arm("survival_only")
de_arm <- arm("de_only")
control <- arm("matched_control")
ridge <- arm("ridge_eligible")
clinical <- arm("clinical")
ridge_sizes <- benchmark_folds$n_model_genes[benchmark_folds$strategy == "ridge_eligible"]
if (any(benchmark$primary_hydra_mismatches != 0L) || nrow(benchmark) != 6L ||
    length(ridge_sizes) != 50L) {
  stop("Nested benchmark summary is incomplete.")
}
benchmark_line <- paste0(
  "The [nested selection benchmark](results/tables/nested_benchmark_summary.csv) ",
  "reused those ", clinical$n_patients, " patients and the same ten-by-five splits. ",
  "The HYDRA arm matched the saved nested-CV gene in every fold. Its concordance ",
  "increment was ", fmt(hydra$mean_delta_c), " (patient-resampled 95% interval ",
  interval(hydra), "). Survival-only selection changed concordance by ",
  fmt(survival_arm$mean_delta_c), " (interval ", interval(survival_arm), ") across ",
  survival_arm$distinct_genes, " genes, and differential-expression-only selection changed it by ",
  fmt(de_arm$mean_delta_c), " (interval ", interval(de_arm), ") across ",
  de_arm$distinct_genes, " genes. Both made the mean three-year Brier score worse. ",
  "An expression-matched control, varying across ", control$distinct_genes,
  " genes, changed concordance by ", fmt(control$mean_delta_c), " (interval ",
  interval(control), "). That interval is above zero and the gain is still several ",
  "times smaller than 0.01, so a small nonzero lower bound is not a useful gain. ",
  "Ridge-penalized Cox regression on the training-fold reproducible genes, averaging ",
  formatC(round(mean(ridge_sizes)), format = "d", big.mark = ","),
  " genes per fold, changed concordance by ", fmt(ridge$mean_delta_c),
  " (interval ", interval(ridge), "). That interval clears 0.01. The mean Brier ",
  "difference was ", fmt(ridge$mean_delta_brier3), " (interval ",
  fmt(ridge$patient_bootstrap_brier_ci_low), " to ",
  fmt(ridge$patient_bootstrap_brier_ci_high),
  "), so the probability-score improvement is not stable under patient resampling. ",
  "Calibration slopes were ", sprintf("%.2f", clinical$mean_calibration_slope),
  " for the clinical model, ", sprintf("%.2f", hydra$mean_calibration_slope),
  " after the HYDRA gene, and ", sprintf("%.2f", ridge$mean_calibration_slope),
  " after the ridge model. The one-gene rules stayed near the clinical model. ",
  "The ridge result means the eligible set still carries held-out ranking information ",
  "when those genes are used together and the penalty is chosen inside the training fold. ",
  "It does not replace the failed one-gene criterion or freeze a signature. ",
  "It was fit only on the TCGA nested splits and was not evaluated on GSE29609 or E-MTAB-1980. ",
  "It was not the prespecified acceptance test, GEO evidence stayed fixed, and the ",
  "patient bootstrap does not refit selection. The side-by-side comparison, including ",
  "ClearCode34, is in [the master figure](results/figures/master_funnel_benchmark.png). ",
  "Repeat-level increments are in [the benchmark figure](results/figures/nested_selection_benchmark.png)."
)
lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
b_start <- which(lines == "<!-- benchmark-results:start -->")
b_end <- which(lines == "<!-- benchmark-results:end -->")
if (length(b_start) != 1L || length(b_end) != 1L || b_end <= b_start) {
  stop("README benchmark-result markers are missing or out of order.")
}
writeLines(c(lines[seq_len(b_start)], benchmark_line, lines[b_end:length(lines)]),
           path, useBytes = TRUE)
message("README nested and benchmark results regenerated from result tables.")
