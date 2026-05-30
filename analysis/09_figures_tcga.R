source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

deg <- read_csv(FILES$tcga_deg, show_col_types = FALSE)
surv <- read_csv(FILES$tcga_survival, show_col_types = FALSE)
candidate_path <- file.path(DIRS$tables, "candidate_gene_evidence_table.csv")
candidates <- if (file.exists(candidate_path)) read_csv(candidate_path, show_col_types = FALSE) else NULL

volcano <- deg |>
  mutate(
    neg_log10_fdr = -log10(padj),
    class = case_when(
      significant & log2FoldChange > 0 ~ "Higher in tumor",
      significant & log2FoldChange < 0 ~ "Lower in tumor",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano, aes(log2FoldChange, neg_log10_fdr, color = class)) +
  geom_point(alpha = 0.45, size = 0.8) +
  geom_vline(xintercept = c(-THRESHOLDS$deg_abs_log2fc, THRESHOLDS$deg_abs_log2fc), linetype = "dashed") +
  geom_hline(yintercept = -log10(THRESHOLDS$deg_fdr), linetype = "dashed") +
  scale_color_manual(values = c("Higher in tumor" = "#B83232", "Lower in tumor" = "#2F6F9F", "Not significant" = "#8A8A8A")) +
  labs(title = "TCGA-KIRC Tumor vs Normal DESeq2", x = "log2 fold change", y = "-log10 FDR", color = NULL) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "tcga_kirc_volcano.png"), p_volcano, width = 7, height = 5, dpi = 300)

if (!is.null(candidates)) {
  discordance <- candidates |>
    filter(!is.na(log_hr), !is.na(fdr), !is.na(tcga_log2fc)) |>
    transmute(
      gene_id = tcga_gene_id,
      symbol = symbol,
      log2FoldChange = tcga_log2fc,
      padj = tcga_padj,
      log_hr = log_hr,
      fdr = fdr,
      abs_log2fc = abs(tcga_log2fc),
      abs_log_hr = abs(log_hr),
      neg_log10_cox_fdr = -log10(fdr),
      prognostic = prognostic,
      reproducible_deg = reproducible_deg
    )
} else {
  discordance <- deg |>
    dplyr::select(gene_id, log2FoldChange, padj) |>
    inner_join(surv, by = "gene_id") |>
    mutate(
      symbol = sub("\\..*$", "", gene_id),
      abs_log2fc = abs(log2FoldChange),
      abs_log_hr = abs(log_hr),
      neg_log10_cox_fdr = -log10(fdr),
      prognostic = fdr < THRESHOLDS$survival_fdr,
      reproducible_deg = FALSE
    )
}

p_discordance <- ggplot(discordance, aes(abs_log2fc, abs_log_hr)) +
  geom_point(aes(size = neg_log10_cox_fdr, fill = prognostic, alpha = reproducible_deg), shape = 21, color = "#222222") +
  geom_text_repel(
    data = discordance |> filter(prognostic, reproducible_deg) |> arrange(fdr) |> slice_head(n = 12),
    aes(label = symbol),
    size = 3,
    max.overlaps = 20
  ) +
  scale_fill_manual(values = c("TRUE" = "#B83232", "FALSE" = "#D9D9D9")) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.22)) +
  labs(
    title = "Differential Expression Magnitude vs Prognostic Effect",
    x = "Absolute tumor-normal log2 fold change",
    y = "Absolute adjusted Cox log(HR)",
    size = "-log10 Cox FDR",
    fill = "Prognostic",
    alpha = "Reproducible"
  ) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "tcga_kirc_discordance.png"), p_discordance, width = 7, height = 5, dpi = 300)

message("Saved TCGA volcano and discordance figures.")
