source("analysis/00_config.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
tab <- function(name) read_csv(file.path(DIRS$tables, paste0(name, ".csv")), show_col_types = FALSE)
coverage <- tab("limitations_plate_coverage")
genes <- tab("limitations_plate_de_all_genes")
candidates <- tab("limitations_plate_de_candidates")
summary <- tab("limitations_plate_de_summary")
priority <- tab("high_confidence_candidate_genes")
stopifnot(nrow(coverage) == 21L, sum(coverage$mixed) == 3L,
          sum(coverage$tumors) == 533L, sum(coverage$normals) == 72L,
          sum(coverage$tumors[coverage$mixed]) == 145L,
          sum(coverage$normals[coverage$mixed]) == 68L,
          nrow(candidates) == 46L, nrow(summary) == 2L)
macro_lines <- readLines("paper/plate_macros.tex")
for (model in c("all_plates", "mixed_plates")) {
  all_rows <- genes |> filter(.data$model == .env$model)
  selected <- candidates |> filter(.data$model == .env$model)
  totals <- summary |> filter(.data$model == .env$model)
  stopifnot(setequal(selected$tcga_gene_id, priority$tcga_gene_id),
            all(selected$same_primary_direction), all(selected$padj < .05),
            setequal(selected$symbol[!selected$significant], c("CRYL1", "RBM47", "KNTC1", "TNFAIP2")),
            all(all_rows$design_rank == all_rows$design_columns),
            all(all_rows$design_columns == all_rows$plates + 1L),
            all(all_rows$n_samples == all_rows$n_tumors + all_rows$n_normals))
  adjusted <- all_rows |> filter(!is.na(padj))
  stopifnot(isTRUE(all.equal(adjusted$padj, p.adjust(adjusted$pvalue, "BH"), tolerance = 1e-10)),
            totals$primary_de_gate_retained == sum(selected$significant & selected$same_primary_direction),
            totals$same_direction == sum(selected$same_primary_direction))
  macro <- if (model == "all_plates") "LimitPlateAllRetained" else "LimitPlateMixedRetained"
  stopifnot(paste0("\\newcommand{\\", macro, "}{", totals$primary_de_gate_retained, "}") %in% macro_lines)
}
message("PASS: plate coverage, full-rank designs, global BH adjustment, fixed shortlist and manuscript counts.")
