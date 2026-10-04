# Verify a memory-only source change before preserving corrected fit results.
source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
old_path <- "results/archive/normalization_20261004/31_frozen_before_memory_cleanup.R"
new_path <- "analysis/31_nested_selection_benchmark.R"
old <- paste(readLines(old_path), collapse = "\n")
expected <- sub("meta <- select_tcga_patient_samples(meta, counts)",
  "meta <- select_tcga_patient_samples(meta, counts)\nrm(se)\ngc()", old, fixed = TRUE)
stopifnot(identical(parse(text = expected), parse(new_path)))
inputs <- c(new_path, "analysis/functions/nested_benchmark.R", "analysis/functions/frozen_normalization.R",
  "analysis/00_config.R", "analysis/functions/patient_samples.R", FILES$tcga_se, FILES$tcga_clinical,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"))
current <- unname(tools::md5sum(inputs))
previous <- current; previous[1] <- unname(tools::md5sum(old_path))
paths <- list.files("data/processed/nested_benchmark_checkpoints", pattern = "^repeat_.*[.]rds$", full.names = TRUE)
backup_dir <- "data/processed/benchmark_memory_cleanup_backups"
dir.create(backup_dir, showWarnings = FALSE)
rows <- list()
for (path in paths) {
  saved <- readRDS(path)
  if (!identical(saved$spec$files, previous)) next
  before <- saved
  backup <- file.path(backup_dir, basename(path))
  if (file.exists(backup)) stop("Backup already exists: ", backup)
  stopifnot(file.copy(path, backup), identical(unname(tools::md5sum(path)), unname(tools::md5sum(backup))))
  saved$spec$files <- current
  write_rds_atomic(saved, path)
  checked <- readRDS(path)
  checked$spec$files <- previous
  stopifnot(identical(checked, before))
  rows[[length(rows) + 1L]] <- tibble(checkpoint = path,
    original_checkpoint_md5 = unname(tools::md5sum(backup)), updated_checkpoint_md5 = unname(tools::md5sum(path)),
    original_source_md5 = previous[1], current_source_md5 = current[1],
    fitted_results_and_rng_unchanged = TRUE,
    scope = "only unused SummarizedExperiment removal and gc added; corrected fitted models retained")
}
stopifnot(length(rows) > 0L)
write_csv_atomic(bind_rows(rows), file.path(DIRS$tables, "benchmark_memory_cleanup_audit.csv"))
message("PASS: ", length(rows), " corrected-source fits preserved; only unused-object removal/gc differ; saved predictions, fold records and RNG unchanged.")
