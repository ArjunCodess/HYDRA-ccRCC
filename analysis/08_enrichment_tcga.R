source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

deg <- read_csv(FILES$tcga_deg, show_col_types = FALSE)
surv <- read_csv(FILES$tcga_survival, show_col_types = FALSE)

gene_ensembl <- surv |>
  filter(model_type == "stage_grade_complete", fdr < THRESHOLDS$strict_survival_fdr) |>
  mutate(ensembl = sub("\\..*$", "", gene_id)) |>
  pull(ensembl) |>
  unique()

if (length(gene_ensembl) < 5) {
  warning("Fewer than 5 prognostic genes passed FDR threshold. Falling back to top 100 survival-ranked genes for exploratory enrichment.")
  gene_ensembl <- surv |>
    filter(model_type == "stage_grade_complete") |>
    arrange(p_value) |>
    slice_head(n = min(100, n())) |>
    mutate(ensembl = sub("\\..*$", "", gene_id)) |>
    pull(ensembl) |>
    unique()
}

universe_ensembl <- deg |>
  filter(!is.na(padj)) |>
  mutate(ensembl = sub("\\..*$", "", gene_id)) |>
  pull(ensembl) |>
  unique()

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ensembl,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) |>
  filter(!is.na(ENTREZID)) |>
  distinct(ENTREZID, .keep_all = TRUE)

universe_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = universe_ensembl,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) |>
  filter(!is.na(ENTREZID)) |>
  distinct(ENTREZID, .keep_all = TRUE)

go_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(universe_map$ENTREZID),
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY")
) |>
  filter(!is.na(GO), ONTOLOGY == "BP") |>
  distinct(ENTREZID, GO)

term_names <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_map$GO),
  keytype = "GOID",
  columns = "TERM"
) |>
  distinct(GOID, .keep_all = TRUE)

gene_set <- intersect(unique(gene_map$ENTREZID), unique(universe_map$ENTREZID))
universe_set <- unique(universe_map$ENTREZID)

M <- length(universe_set)
N <- length(gene_set)

ego_df <- go_map |>
  group_by(GO) |>
  summarise(
    term_size = n_distinct(ENTREZID),
    overlap = sum(unique(ENTREZID) %in% gene_set),
    genes = paste(sort(unique(ENTREZID[ENTREZID %in% gene_set])), collapse = ";"),
    .groups = "drop"
  ) |>
  filter(overlap > 0, term_size >= 10, term_size <= 500) |>
  mutate(
    pvalue = phyper(overlap - 1, term_size, M - term_size, N, lower.tail = FALSE),
    p.adjust = stats::p.adjust(pvalue, method = "BH"),
    gene_ratio = overlap / N,
    background_ratio = term_size / M
  ) |>
  left_join(term_names, by = c("GO" = "GOID")) |>
  arrange(p.adjust, pvalue) |>
  dplyr::select(ID = GO, Description = TERM, term_size, overlap, gene_ratio, background_ratio, pvalue, p.adjust, genes)

write_csv_atomic(ego_df, FILES$tcga_enrichment)

message("GO BP enrichment complete. Terms saved: ", nrow(ego_df))
