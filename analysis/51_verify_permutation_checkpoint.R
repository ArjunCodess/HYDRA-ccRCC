# One-time validation-only cache repair; no statistical computation is changed.
source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
old_path <- "results/archive/normalization_20261004/33_checkpointed_source.R"
new_path <- "analysis/33_survival_permutation_null.R"
old_source <- paste(readLines(old_path), collapse = "\n")
expected_source <- gsub("identical(result$gene_id, original$gene_id)",
  "identical(as.character(result$gene_id), as.character(original$gene_id))", old_source, fixed = TRUE)
expected_source <- gsub("all.equal(result[[column]], original[[column]], tolerance = 1e-10)",
  "all.equal(as.numeric(result[[column]]), as.numeric(original[[column]]), tolerance = 1e-10)",
  expected_source, fixed = TRUE)
stopifnot(identical(parse(text = expected_source), parse(new_path)))
folds <- sprintf("data/processed/nested_cv_fold_checkpoints/repeat_01_fold_%d.rds", 1:5)
inputs <- c(new_path, "analysis/00_config.R", "analysis/functions/io.R",
  "analysis/functions/tcga_metadata.R", "analysis/functions/patient_samples.R",
  "environment/package_versions.csv", FILES$tcga_se, FILES$tcga_clinical, folds)
current <- unname(tools::md5sum(inputs))
path <- "data/processed/review_permutation_iterations.rds"
saved <- readRDS(path)
stopifnot(saved$signature$permutations == 200L,
          identical(saved$signature$r_version, R.version.string),
          identical(saved$signature$files[-1], current[-1]),
          identical(saved$signature$files[1], unname(tools::md5sum(old_path))),
          length(saved$rows) == 200L)
result <- bind_rows(saved$rows)
original <- read_csv("results/archive/normalization_20261004/survival_permutation_null.csv", show_col_types = FALSE)
stopifnot(identical(result$simulation, 1:200),
          identical(as.character(result$gene_id), as.character(original$gene_id)))
for (column in c("main_fdr", "main_beta")) stopifnot(isTRUE(all.equal(
  as.numeric(result[[column]]), as.numeric(original[[column]]), tolerance = 1e-10)))
old_md5 <- unname(tools::md5sum(path))
saved$signature$files <- current
write_rds_atomic(saved, path)
restored <- readRDS(path)
stopifnot(identical(restored$rows, saved$rows), identical(restored$rng_after, saved$rng_after))
write_csv_atomic(tibble(iterations = 200L, original_gene_selections_reproduced = TRUE,
  statistical_source_unchanged = TRUE, old_source_md5 = unname(tools::md5sum(old_path)),
  corrected_source_md5 = current[1], original_checkpoint_md5 = old_md5,
  repaired_checkpoint_md5 = unname(tools::md5sum(path)),
  scope = "validation-only type coercion for all-missing CSV columns; rows and RNG unchanged"),
  file.path(DIRS$tables, "permutation_checkpoint_type_audit.csv"))
message("PASS: only comparison types changed; all 200 selections and coefficients match; checkpoint rows/RNG preserved.")
