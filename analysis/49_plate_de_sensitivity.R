# Technical-plate proxy sensitivity. Neither model replaces the primary DE gate.
source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/patient_samples.R")
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(dplyr); library(readr)})
set.seed(RESAMPLING$seed + 49L)
se <- read_required_rds(FILES$tcga_se)
counts <- assay(se, "unstranded")
meta <- as.data.frame(colData(se)) |> tibble::rownames_to_column("sample_barcode") |>
  mutate(condition = factor(shortLetterCode, levels = c("NT", "TP")))
meta <- select_tcga_patient_samples(meta, counts) |>
  mutate(plate = factor(substr(sample_barcode, 22L, 25L)))
rm(se); gc()
stopifnot(nrow(meta) == 605L, all(nchar(meta$sample_barcode) == 28L))
plate_coverage <- meta |> group_by(plate) |> summarise(tumors = sum(condition == "TP"),
  normals = sum(condition == "NT"), mixed = tumors > 0 & normals > 0, .groups = "drop")
write_csv_atomic(plate_coverage, file.path(DIRS$tables, "limitations_plate_coverage.csv"))
priority <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), show_col_types = FALSE)
primary <- read_csv(FILES$tcga_deg, show_col_types = FALSE)
signature <- list(seed = RESAMPLING$seed + 49L,
  inputs = unname(tools::md5sum(c("analysis/49_plate_de_sensitivity.R", "analysis/00_config.R",
    "analysis/functions/patient_samples.R", FILES$tcga_se, FILES$tcga_deg,
    file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), "environment/package_versions.csv"))))
results <- list(); summaries <- list()
for (model in c("all_plates", "mixed_plates")) {
  metadata <- if (model == "all_plates") meta else meta |> filter(as.character(plate) %in% as.character(plate_coverage$plate[plate_coverage$mixed]))
  metadata <- as.data.frame(metadata); metadata$plate <- droplevels(metadata$plate)
  rownames(metadata) <- metadata$sample_barcode
  design <- model.matrix(~plate + condition, metadata)
  stopifnot(qr(design)$rank == ncol(design), sum(metadata$condition == "NT") > 0,
            sum(metadata$condition == "TP") > 0)
  raw <- counts[, metadata$sample_barcode, drop = FALSE]
  keep <- rowSums(raw >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
  checkpoint <- file.path(DIRS$processed, paste0("limitations_plate_", model, ".rds"))
  cached <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
  if (!is.null(cached) && identical(cached$signature, signature)) {
    result <- cached$result; fallback <- cached$fallback
    message("Reused verified plate sensitivity: ", model)
  } else {
    message("Fitting plate DE sensitivity: ", model, "; ", nrow(metadata), " samples, ", ncol(design), " parameters")
    dds <- DESeqDataSetFromMatrix(raw[keep, , drop = FALSE], metadata, design = ~plate + condition)
    fitted <- tryCatch(DESeq(dds, quiet = TRUE), error = identity)
    fallback <- inherits(fitted, "error")
    if (fallback) {
      if (!identical(conditionMessage(fitted), "default method not implemented for type 'expression'") ||
          !identical(deparse(conditionCall(fitted)), "is.finite(partial)")) stop(fitted)
      message("Pinned outlier-replacement error; retrying without replacement.")
      fitted <- DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    }
    result <- as.data.frame(DESeq2::results(fitted, contrast = c("condition", "TP", "NT"))) |>
      tibble::rownames_to_column("gene_id")
    write_rds_atomic(list(signature = signature, result = result, fallback = fallback), checkpoint)
    rm(dds, fitted); gc()
  }
  result <- result |> mutate(model, n_samples = nrow(metadata), n_tumors = sum(metadata$condition == "TP"),
    n_normals = sum(metadata$condition == "NT"), plates = nlevels(metadata$plate),
    design_rank = qr(design)$rank, design_columns = ncol(design), outlier_replacement_fallback = fallback,
    significant = !is.na(padj) & padj < THRESHOLDS$deg_fdr & abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc)
  results[[model]] <- result
  candidate <- priority |> select(symbol, tcga_gene_id) |>
    left_join(result, by = c("tcga_gene_id" = "gene_id")) |>
    left_join(primary |> select(gene_id, primary_log2fc = log2FoldChange), by = c("tcga_gene_id" = "gene_id")) |>
    mutate(same_primary_direction = sign(log2FoldChange) == sign(primary_log2fc))
  stopifnot(nrow(candidate) == 23L, all(is.finite(candidate$log2FoldChange)))
  summaries[[model]] <- candidate
}
write_csv_atomic(bind_rows(results), file.path(DIRS$tables, "limitations_plate_de_all_genes.csv"))
candidate_results <- bind_rows(summaries)
write_csv_atomic(candidate_results, file.path(DIRS$tables, "limitations_plate_de_candidates.csv"))
summary <- candidate_results |> group_by(model) |>
  summarise(n_samples = first(n_samples), n_tumors = first(n_tumors), n_normals = first(n_normals),
    plates = first(plates), same_direction = sum(same_primary_direction),
    fdr_supported = sum(!is.na(padj) & padj < .05), primary_de_gate_retained = sum(significant & same_primary_direction),
    tested_genes = nrow(results[[first(model)]]), .groups = "drop")
write_csv_atomic(summary, file.path(DIRS$tables, "limitations_plate_de_summary.csv"))
writeLines(capture.output(sessionInfo()), "environment/plate_sensitivity_sessionInfo.txt")
macro <- function(name, value) paste0("\\newcommand{\\", name, "}{", value, "}")
value <- function(model, column) summary[[column]][summary$model == model]
writeLines(c("% Generated by analysis/49_plate_de_sensitivity.R.",
  macro("LimitPlateAllRetained", value("all_plates", "primary_de_gate_retained")),
  macro("LimitPlateMixedRetained", value("mixed_plates", "primary_de_gate_retained")),
  macro("LimitMixedTumors", value("mixed_plates", "n_tumors")),
  macro("LimitMixedNormals", value("mixed_plates", "n_normals"))), "paper/plate_macros.tex")
message("Plate sensitivities complete; no candidate gates redefined.")
