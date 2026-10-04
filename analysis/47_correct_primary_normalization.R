# Audited held-out transformation correction. Preserve training selections and
# reproduce every original prediction before migrating any checkpoint.
archive <- "results/archive/normalization_20261004"
old_script <- file.path(archive, "22_nested_cv.R")
stopifnot(file.exists(old_script))
old_lines <- readLines(old_script)
new_lines <- readLines("analysis/22_nested_cv.R")
# Only the normalization call, helper import and signature dependency may differ.
canonical <- function(lines) {
  lines <- lines[!grepl('source\\("analysis/functions/frozen_normalization.R"\\)|^ +"analysis/functions/frozen_normalization.R",', lines)]
  lines <- gsub("frozen_size_factors(test_raw, geo_mean)",
                "DESeq2::estimateSizeFactorsForMatrix(test_raw, geoMeans = geo_mean)", lines, fixed = TRUE)
  parse(text = lines, keep.source = FALSE)
}
stopifnot(identical(canonical(old_lines), canonical(new_lines)))
# Load setup and functions, never the outer loop or writers.
expressions <- parse("analysis/22_nested_cv.R")
for (expression in expressions) {
  if (identical(expression[[1]], as.name("<-")) && identical(expression[[2]], as.name("predictions"))) break
  eval(expression, envir = .GlobalEnv)
}
old_inputs <- c(old_script, "analysis/00_config.R", "analysis/functions/patient_samples.R",
  FILES$tcga_se, FILES$tcga_clinical,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"))
old_signature <- unname(tools::md5sum(old_inputs))
new_inputs <- c("analysis/22_nested_cv.R", "analysis/00_config.R",
  "analysis/functions/frozen_normalization.R", old_inputs[3:7])
new_signature <- unname(tools::md5sum(new_inputs))
records <- list(); migrated <- list()
for (repeat_id in 1:10) {
  set.seed(RESAMPLING$seed + 22L + repeat_id)
  fold_id <- folds(patients$os_event, 5L)
  for (fold in 1:5) {
    name <- sprintf("repeat_%02d_fold_%d.rds", repeat_id, fold)
    source_path <- file.path(archive, "nested_cv_fold_checkpoints", name)
    saved <- readRDS(source_path)
    stopifnot(identical(saved$signature, old_signature))
    train_idx <- which(fold_id != fold); test_idx <- which(fold_id == fold)
    normals <- all_normal |> filter(!patient_barcode %in% patients$patient_barcode[test_idx])
    training_samples <- c(patients$sample_barcode[train_idx], normals$sample_barcode)
    raw <- counts[, training_samples, drop = FALSE]
    raw <- raw[rowSums(raw >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples, , drop = FALSE]
    reference <- exp(rowMeans(log(raw)))
    train_sf <- DESeq2::estimateSizeFactorsForMatrix(raw)
    corrected_sf <- frozen_size_factors(tumor_counts[rownames(raw), test_idx, drop = FALSE], reference)
    batch_sf <- DESeq2::estimateSizeFactorsForMatrix(tumor_counts[rownames(raw), test_idx, drop = FALSE], geoMeans = reference)
    gene <- saved$summary$selected_gene
    ids <- if (!is.null(saved$cache)) rownames(saved$cache$fd$expr) else gene[!is.na(gene)]
    normalized <- function(test_sf) {
      sf <- c(setNames(train_sf, training_samples), setNames(test_sf, patients$sample_barcode[test_idx]))
      log2(sweep(counts[ids, patients$sample_barcode, drop = FALSE], 2L, sf[patients$sample_barcode], "/") + 1)
    }
    old_expr <- normalized(batch_sf); corrected_expr <- normalized(corrected_sf)
    if (!is.null(saved$cache)) {
      stopifnot(max(abs(old_expr - saved$cache$fd$expr)) < 1e-10)
    }
    original <- predict_fold(list(expr = old_expr), train_idx, test_idx, gene, patients$os_time, patients$os_event) |>
      mutate(repeat_id = repeat_id, fold = fold)
    columns <- c("clinical_lp", "gene_lp", "clinical_risk3", "gene_risk3")
    stopifnot(identical(original$patient_barcode, saved$prediction$patient_barcode))
    error <- max(abs(as.matrix(original[, columns]) - as.matrix(saved$prediction[, columns])))
    stopifnot(is.finite(error), error < 1e-9)
    corrected <- predict_fold(list(expr = corrected_expr), train_idx, test_idx, gene, patients$os_time, patients$os_event) |>
      mutate(repeat_id = repeat_id, fold = fold)
    stopifnot(identical(original$clinical_lp, corrected$clinical_lp),
              max(abs(old_expr[, train_idx, drop = FALSE] - corrected_expr[, train_idx, drop = FALSE])) == 0)
    records[[length(records) + 1L]] <- tibble(repeat_id, fold, source_path,
      original_md5 = unname(tools::md5sum(source_path)), original_prediction_max_error = error,
      training_expression_max_change = 0, test_factor_recentering = exp(mean(log(corrected_sf))),
      max_gene_lp_change = max(abs(corrected$gene_lp-original$gene_lp)),
      selected_gene = gene, correction = "test transform only; unchanged DE and selection; Cox refit")
    saved$prediction <- corrected; saved$signature <- new_signature
    if (!is.null(saved$cache)) saved$cache$fd$expr <- corrected_expr
    saved$normalization_correction <- records[[length(records)]]
    migrated[[name]] <- saved
    message("Verified original and corrected primary repeat ", repeat_id, " fold ", fold)
  }
}
# No mutation until all 50 original predictions and signatures have passed.
stopifnot(length(migrated) == 50L)
for (name in names(migrated)) write_rds_atomic(migrated[[name]], file.path(DIRS$processed, "nested_cv_fold_checkpoints", name))
write_csv_atomic(bind_rows(records), file.path(DIRS$tables, "normalization_primary_correction_audit.csv"))
writeLines(capture.output(sessionInfo()), "environment/normalization_correction_sessionInfo.txt")
message("PASS: 50 original predictions reproduced; corrected test transforms; training and selection unchanged.")
