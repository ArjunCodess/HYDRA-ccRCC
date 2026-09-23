source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages(library(readr))

nested <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"), show_col_types = FALSE)
funnel <- read_csv(file.path(DIRS$tables, "funnel_external_paired_bootstrap.csv"),
                   show_col_types = FALSE)
stopifnot(nrow(nested) == 1L, nested$repeats == 10L, nested$folds == 5L,
          nrow(funnel) == 2L)

criteria <- data.frame(
  criterion = c("selection_aware_cv", "external_funnel_vs_survival_only",
                "untouched_external_patients", "tumor_single_cell_evidence",
                "independent_randomized_treatment_cohort"),
  status = c(ifelse(nested$cv_acceptance, "pass", "fail"),
             ifelse(all(funnel$ci_low > 0), "pass", "fail"),
             "unavailable", "unavailable", "unavailable"),
  evidence = c("nested_cv_summary.csv; nested_cv_patient_bootstrap.csv",
               "funnel_external_paired_bootstrap.csv", "previously inspected cohorts only",
               "HPA normal tissue only", "CheckMate 025 is post hoc and already inspected"),
  rule = c("mean delta C >= 0.01; patient bootstrap lower 95% > 0; mean delta Brier <= 0",
           "complete minus survival-only directional rate lower paired 95% > 0 in both inspected cohorts",
           "new patients not yet inspected", "ccRCC tumor-cell-resolved expression",
           "prespecified independent treatment interaction"),
  stringsAsFactors = FALSE
)
write_csv_atomic(criteria, file.path(DIRS$tables, "acceptance_criteria.csv"))
message("Scientific acceptance: ", paste(criteria$status[1:2], collapse = ", "), ".")
