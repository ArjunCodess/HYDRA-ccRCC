source("analysis/00_config.R")
suppressPackageStartupMessages(library(readr))

nested <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"),
                   show_col_types = FALSE)
null <- read_csv(file.path(DIRS$tables, "nested_cv_clinical_null_summary.csv"),
                 show_col_types = FALSE)
if (nrow(nested) != 1L || nested$repeats != 10L || nested$folds != 5L ||
    nrow(null) != 1L || null$simulations != 200L) {
  stop("Final nested CV and null summaries are required before updating README.")
}

line <- paste0(
  "The [selection-aware nested CV](results/tables/nested_cv_summary.csv) used ",
  nested$n_patients, " patients across ten repeats of five folds. It selected no gene in ",
  nested$no_gene_folds, " folds; ", nested$outlier_replacement_fallbacks,
  " folds required a DESeq2 no-replacement retry. The mean selected-gene minus ",
  "clinical concordance was ", sprintf("%+.4f", nested$mean_delta_c),
  " (patient-resampled 95% interval ",
  sprintf("%+.4f", nested$patient_bootstrap_ci_low), " to ",
  sprintf("%+.4f", nested$patient_bootstrap_ci_high),
  "), and the three-year Brier-score difference was ",
  sprintf("%+.4f", nested$mean_delta_brier3), ". The prespecified CV criterion **",
  ifelse(nested$cv_acceptance, "passed", "failed"),
  "**. The [200 clinical-only null simulations](results/tables/nested_cv_clinical_null_summary.csv) ",
  "selected a gene in ", null$gene_selected,
  " single-fold runs. The patient bootstrap keeps the fitted fold models fixed, ",
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
message("README nested results regenerated from committed table schemas.")
