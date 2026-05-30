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
    filter(!is.na(main_log_hr), !is.na(main_fdr), !is.na(tcga_log2fc)) |>
    transmute(
      gene_id = tcga_gene_id,
      symbol = symbol,
      log2FoldChange = tcga_log2fc,
      padj = tcga_padj,
      log_hr = main_log_hr,
      fdr = main_fdr,
      abs_log2fc = abs(tcga_log2fc),
      abs_log_hr = abs(main_log_hr),
      neg_log10_cox_fdr = -log10(main_fdr),
      prognostic = prognostic,
      reproducible_deg = reproducible_deg,
      strict_candidate = strict_candidate,
      high_confidence_candidate = high_confidence_candidate,
      pathway_class = pathway_class
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
      prognostic = fdr < THRESHOLDS$strict_survival_fdr,
      reproducible_deg = FALSE,
      strict_candidate = FALSE,
      high_confidence_candidate = FALSE,
      pathway_class = "Unclassified"
  )
}

p_directional <- ggplot(discordance, aes(log2FoldChange, log_hr)) +
  geom_hline(yintercept = 0, color = "#5A5A5A", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#5A5A5A", linewidth = 0.35) +
  geom_point(aes(size = neg_log10_cox_fdr, fill = pathway_class, alpha = reproducible_deg, stroke = high_confidence_candidate), shape = 21, color = "#222222") +
  geom_text_repel(
    data = discordance |> filter(symbol %in% c("KL", "ACADM", "CRYL1", "ACAT1", "DDC", "TCIRG1", "HHLA2", "C1QTNF6")),
    aes(label = symbol),
    size = 3,
    max.overlaps = 20
  ) +
  scale_fill_manual(values = c(
    "Hypoxia" = "#7B3294",
    "Angiogenesis" = "#008837",
    "ECM/EMT" = "#A6611A",
    "Metabolism" = "#0571B0",
    "Immune" = "#C51B7D",
    "Hypoxia;Metabolism" = "#5E3C99",
    "ECM/EMT;Immune" = "#B35806",
    "Immune;Metabolism" = "#E66101",
    "Unclassified" = "#C9C9C9"
  ), na.value = "#C9C9C9") +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.22)) +
  scale_discrete_manual(aesthetics = "stroke", values = c("TRUE" = 1.1, "FALSE" = 0.25), guide = "none") +
  labs(
    title = "Direction of Tumor-Normal Change vs Prognostic Effect",
    x = "TCGA tumor-normal log2 fold change",
    y = "Stage/grade-adjusted Cox log(HR)",
    size = "-log10 main Cox FDR",
    fill = "Pathway class",
    alpha = "Reproducible"
  ) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "tcga_kirc_directional_discordance.png"), p_directional, width = 7, height = 5, dpi = 300)

p_discordance <- ggplot(discordance, aes(abs_log2fc, abs_log_hr)) +
  geom_point(aes(size = neg_log10_cox_fdr, fill = pathway_class, alpha = reproducible_deg, stroke = high_confidence_candidate), shape = 21, color = "#222222") +
  geom_text_repel(
    data = discordance |> filter(symbol %in% c("KL", "ACADM", "CRYL1", "ACAT1", "DDC", "TCIRG1", "HHLA2", "C1QTNF6")),
    aes(label = symbol),
    size = 3,
    max.overlaps = 20
  ) +
  scale_fill_manual(values = c(
    "Hypoxia" = "#7B3294",
    "Angiogenesis" = "#008837",
    "ECM/EMT" = "#A6611A",
    "Metabolism" = "#0571B0",
    "Immune" = "#C51B7D",
    "Hypoxia;Metabolism" = "#5E3C99",
    "ECM/EMT;Immune" = "#B35806",
    "Immune;Metabolism" = "#E66101",
    "Unclassified" = "#C9C9C9"
  ), na.value = "#C9C9C9") +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.22)) +
  scale_discrete_manual(aesthetics = "stroke", values = c("TRUE" = 1.1, "FALSE" = 0.25), guide = "none") +
  labs(
    title = "Reproducible Expression Magnitude vs Strict Prognostic Effect",
    x = "Absolute tumor-normal log2 fold change",
    y = "Absolute stage/grade-adjusted Cox log(HR)",
    size = "-log10 main Cox FDR",
    fill = "Pathway class",
    alpha = "Reproducible"
  ) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "tcga_kirc_discordance.png"), p_discordance, width = 7, height = 5, dpi = 300)

candidate_priority_path <- file.path(DIRS$tables, "manuscript_candidate_prioritization.csv")
if (file.exists(candidate_priority_path)) {
  candidate_priority <- read_csv(candidate_priority_path, show_col_types = FALSE) |>
    mutate(
      manual_tier = factor(
        manual_tier,
        levels = c("lead", "supporting", "supporting risk", "interpret cautiously", "composition flag", "do not highlight")
      ),
      symbol_label = paste0(symbol, " (", manual_tier, ")"),
      symbol_label = factor(symbol_label, levels = rev(symbol_label)),
      log_hr_ci_low = log(main_hr_ci_low),
      log_hr_ci_high = log(main_hr_ci_high)
    )

  p_forest <- ggplot(candidate_priority, aes(main_log_hr, symbol_label, color = manual_tier)) +
    geom_vline(xintercept = 0, color = "#555555", linewidth = 0.4) +
    geom_errorbar(aes(xmin = log_hr_ci_low, xmax = log_hr_ci_high), orientation = "y", width = 0, linewidth = 0.45) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c(
      "lead" = "#1B7837",
      "supporting" = "#0571B0",
      "supporting risk" = "#B35806",
      "interpret cautiously" = "#7B3294",
      "composition flag" = "#777777",
      "do not highlight" = "#B83232"
    ), drop = FALSE) +
    labs(
      title = "High-Confidence Candidates Require Unequal Manuscript Weight",
      x = "Stage/grade-adjusted Cox log(HR) per SD expression",
      y = NULL,
      color = "Manuscript tier"
    ) +
    theme_hydra() +
    theme(legend.position = "bottom")

  ggsave(file.path(DIRS$figures, "candidate_forest_plot.png"), p_forest, width = 7, height = 6.4, dpi = 300)
}

message("Saved TCGA volcano, discordance, and candidate figures.")
