source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages(library(readr))

null <- read_csv(file.path(DIRS$tables, "nested_cv_clinical_null.csv"),
                 show_col_types = FALSE)
if (nrow(null) != 200L || anyDuplicated(null$simulation) ||
    any(!is.finite(null$delta_c))) {
  stop("Clinical-only null simulations are incomplete.")
}
summary <- data.frame(
  simulations = nrow(null),
  gene_selected = sum(!is.na(null$selected_gene)),
  selection_fraction = mean(!is.na(null$selected_gene)),
  mean_fold_delta_c = mean(null$delta_c),
  median_fold_delta_c = median(null$delta_c),
  q95_fold_delta_c = unname(quantile(null$delta_c, 0.95))
)
write_csv_atomic(summary, file.path(DIRS$tables, "nested_cv_clinical_null_summary.csv"))
message("Summarized ", nrow(null), " clinical-only null simulations.")
