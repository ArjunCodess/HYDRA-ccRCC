source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

candidates <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
                       show_col_types = FALSE) |>
  select(symbol, tcga_gene_id, primary_log2fc = tcga_log2fc)
paired <- read_csv(file.path(DIRS$tables, "tcga_kirc_paired_deseq2_tumor_vs_normal.csv"),
                   show_col_types = FALSE) |>
  select(gene_id, paired_log2fc = log2FoldChange, paired_padj = padj,
         paired_significant = significant)
out <- candidates |> left_join(paired, by = c("tcga_gene_id" = "gene_id")) |>
  mutate(same_direction = sign(primary_log2fc) == sign(paired_log2fc))
stopifnot(nrow(out) == nrow(candidates), !anyDuplicated(out$tcga_gene_id))
write_csv_atomic(out, file.path(DIRS$tables, "candidate_paired_de_sensitivity.csv"))
message("Paired DE retained direction for ", sum(out$same_direction, na.rm = TRUE),
        " candidates; ", sum(out$paired_significant, na.rm = TRUE),
        " met the paired DE threshold.")
