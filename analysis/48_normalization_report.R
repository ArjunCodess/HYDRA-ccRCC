source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
archive <- "results/archive/normalization_20261004"
old_primary <- read_csv(file.path(archive, "nested_cv_summary.csv"), show_col_types = FALSE)
new_primary <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"), show_col_types = FALSE)
old_benchmark <- read_csv(file.path(archive, "nested_benchmark_summary.csv"), show_col_types = FALSE)
new_benchmark <- read_csv(file.path(DIRS$tables, "nested_benchmark_summary.csv"), show_col_types = FALSE)
# No report before all benchmark checkpoints have the corrected-source signature.
checkpoint_inputs <- c("analysis/31_nested_selection_benchmark.R", "analysis/functions/nested_benchmark.R",
  "analysis/functions/frozen_normalization.R", "analysis/00_config.R", "analysis/functions/patient_samples.R",
  FILES$tcga_se, FILES$tcga_clinical, file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"))
expected <- list(files = unname(tools::md5sum(checkpoint_inputs)),
  screen_n = as.integer(Sys.getenv("HYDRA_SURVIVAL_SCREEN", "2000")), seed = RESAMPLING$seed,
  ridge = "glmnet-5.0-efron-alpha0-lambda.min")
paths <- list.files(file.path(DIRS$processed, "nested_benchmark_checkpoints"), pattern = "^repeat_.*[.]rds$", full.names = TRUE)
stopifnot(length(paths) == 50L)
for (path in paths) stopifnot(identical(readRDS(path)$spec, expected))
combine <- function(old, current, evaluation) {
  if (!"strategy" %in% names(old)) old$strategy <- current$strategy <- "hydra"
  columns <- c("strategy", "mean_delta_c", "mean_delta_brier3", "patient_bootstrap_ci_low", "patient_bootstrap_ci_high")
  inner_join(old[, columns], current[, columns], by = "strategy", suffix = c("_before", "_corrected")) |>
    mutate(evaluation, delta_c_change = mean_delta_c_corrected - mean_delta_c_before,
           uncertainty_scope = "conditional bootstrap of saved predictions; excludes training and selection uncertainty")
}
comparison <- bind_rows(combine(old_primary, new_primary, "primary"), combine(old_benchmark, new_benchmark, "benchmark"))
stopifnot(nrow(comparison) == 7L,
          abs(new_primary$mean_delta_c - new_benchmark$mean_delta_c[new_benchmark$strategy == "hydra"]) < 1e-10)
write_csv_atomic(comparison, file.path(DIRS$tables, "normalization_prediction_comparison.csv"))
audit <- read_csv(file.path(DIRS$tables, "normalization_primary_correction_audit.csv"), show_col_types = FALSE)
benchmark_fits <- bind_rows(lapply(paths, function(path) {
  saved <- readRDS(path)
  tibble(checkpoint = path, md5 = unname(tools::md5sum(path)),
    repeat_id = saved$fold_rows$repeat_id[1], fold = saved$fold_rows$fold[1],
    elapsed_seconds = saved$fold_rows$elapsed_seconds[1],
    ridge_fallback_to_clinical = saved$fold_rows$fallback_to_clinical[saved$fold_rows$strategy == "ridge_eligible"],
    outlier_replacement_fallback = saved$fold_rows$outlier_replacement_fallback[1],
    status = "fresh fit under frozen-reference normalization; memory-only source equivalence recorded separately")
}))
stopifnot(!anyDuplicated(paste(benchmark_fits$repeat_id, benchmark_fits$fold)),
          all(table(benchmark_fits$repeat_id) == 5L))
write_csv_atomic(benchmark_fits, file.path(DIRS$tables, "normalization_benchmark_refit_audit.csv"))
format_signed <- function(x) sprintf("%+.4f", x)
labels <- new_benchmark |> filter(strategy != "clinical") |>
  transmute(strategy, delta_c = format_signed(mean_delta_c),
            interval = paste0("[", format_signed(patient_bootstrap_ci_low), ", ", format_signed(patient_bootstrap_ci_high), "]"))
label_rows <- apply(as.data.frame(labels), 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
svg <- "paper/figures/hydra_evidence_overview.svg"
writeLines(c("# Author-maintained overview", "",
  "The author's SVG remains Figure 1 and the first README image. Stage 45 exports it without changing its bytes, layout, or labels. The author-updated artwork predates the 2026-10-04 normalization correction, so its prediction labels and marks need the edits below. Other funnel counts, external direction counts, and composition counts remain current.", "",
  "## Prediction labels to update", "",
  "Replace the one-gene and affected comparator point labels and move their plotted marks to the corrected positions. Use the benchmark interval on the HYDRA whisker; the primary interval uses a different bootstrap seed. Keep the original design.", "",
  "| Benchmark arm | Corrected concordance increment | Conditional 95% interval |",
  "| --- | --- | --- |", label_rows, "",
  "ClearCode34 uses per-sample CPM and is unaffected by this correction. Keep its existing point. These are research-procedure comparisons, not a validated panel.", "",
  sprintf("The corrected primary HYDRA estimate is %s [%s, %s]; the one-gene criterion still fails.",
    format_signed(new_primary$mean_delta_c), format_signed(new_primary$patient_bootstrap_ci_low), format_signed(new_primary$patient_bootstrap_ci_high)), "",
  "## Wording to clarify", "",
  "Label the external direction counts as unadjusted; retain 21/22 and 5/21. Harmonized clinical adjustment is a separate assessment reported in the paper.", "",
  "Replace 'Broader reproducible gene space retains held-out ranking signal.' with 'The joint ridge predictor improves held-out ranking within TCGA.' Ridge estimates clinical and gene coefficients together, so its contrast does not isolate the contribution of individual genes.", "",
  "## Export and source evidence", "",
  "Run `python analysis/45_evidence_overview.py` after an author edit, then rebuild the paper. The export uses Chrome/Edge/Chromium and the SVG viewBox dimensions. Temporary print HTML and the isolated browser profile stay in ignored `paper/qa_render/`.", "",
  "Values above are generated by `analysis/48_normalization_report.R` from `results/tables/nested_benchmark_summary.csv` and `nested_cv_summary.csv`. `normalization_prediction_comparison.csv` preserves before/after estimates, `normalization_primary_correction_audit.csv` verifies unchanged training and original predictions, and `normalization_benchmark_refit_audit.csv` identifies all 50 fresh benchmark fits. Original outputs are in `results/archive/normalization_20261004`.", "",
  sprintf("Predictions from all %d original primary folds were reproduced before correction; maximum linear-predictor change was %.6f. The artwork is preserved and its prediction values are identified as pre-correction in the manuscript caption and README.", nrow(audit), max(audit$max_gene_lp_change))),
  "docs/HYDRA_OVERVIEW_FIGURE_EDIT_GUIDE.md")
message("PASS: 50 fresh benchmark fits; before/after prediction evidence; author SVG edits documented, not applied.")
