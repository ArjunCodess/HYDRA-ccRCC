source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

tcga <- read_csv(FILES$tcga_deg, show_col_types = FALSE) |>
  mutate(ensembl = sub("\\..*$", "", gene_id))

tcga_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(tcga$ensembl),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) |>
  filter(!is.na(SYMBOL)) |>
  distinct(ENSEMBL, .keep_all = TRUE)

tcga_symbol <- tcga |>
  left_join(tcga_map, by = c("ensembl" = "ENSEMBL")) |>
  filter(!is.na(SYMBOL)) |>
  group_by(SYMBOL) |>
  arrange(padj, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    symbol = SYMBOL,
    tcga_gene_id = gene_id,
    tcga_log2fc = log2FoldChange,
    tcga_padj = padj,
    tcga_significant = !is.na(padj) & padj < THRESHOLDS$deg_fdr & abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc,
    tcga_direction = sign(log2FoldChange)
  )

read_geo <- function(accession) {
  prefix <- tolower(accession)
  geo <- read_csv(file.path(DIRS$tables, paste0(prefix, "_limma_tumor_vs_normal.csv")), show_col_types = FALSE)
  out <- geo |>
    transmute(
      symbol = symbol,
      log2fc = log2FoldChange,
      pvalue = pvalue,
      padj = padj,
      direction = sign(log2fc),
      nominal = pvalue < 0.05
    )
  names(out) <- c("symbol", paste0(prefix, c("_log2fc", "_pvalue", "_padj", "_direction", "_nominal")))
  out
}

gse40435 <- read_geo("GSE40435")
gse53757 <- read_geo("GSE53757")

repro <- tcga_symbol |>
  left_join(gse40435, by = "symbol") |>
  left_join(gse53757, by = "symbol") |>
  mutate(
    same_direction_gse40435 = !is.na(gse40435_direction) & gse40435_direction == tcga_direction,
    same_direction_gse53757 = !is.na(gse53757_direction) & gse53757_direction == tcga_direction,
    same_direction_count = rowSums(cbind(same_direction_gse40435, same_direction_gse53757), na.rm = TRUE),
    nominal_support_count = rowSums(cbind(
      same_direction_gse40435 & gse40435_nominal,
      same_direction_gse53757 & gse53757_nominal
    ), na.rm = TRUE),
    reproducible_deg = tcga_significant & same_direction_count >= 2 & nominal_support_count >= 1
  ) |>
  arrange(desc(reproducible_deg), tcga_padj)

write_csv_atomic(repro, file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"))

summary <- tibble(
  metric = c("tcga_symbol_genes", "tcga_significant", "same_direction_both_geo", "reproducible_deg"),
  value = c(
    nrow(tcga_symbol),
    sum(tcga_symbol$tcga_significant, na.rm = TRUE),
    sum(repro$same_direction_count >= 2, na.rm = TRUE),
    sum(repro$reproducible_deg, na.rm = TRUE)
  )
)
write_csv_atomic(summary, file.path(DIRS$tables, "reproducibility_summary.csv"))

message("Reproducible DEG count: ", sum(repro$reproducible_deg, na.rm = TRUE))
