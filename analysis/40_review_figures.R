# Main-paper figures for the computational-genomics review revision.
# Run from repository root: Rscript analysis/40_review_figures.R
# Inputs are regenerated study result tables, never hand-entered effect estimates.
source("analysis/00_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})
tab <- function(name) read_csv(file.path(DIRS$tables, name), show_col_types = FALSE)
sources <- c("evidence_funnel.csv", "external_survival_emtab1980_summary.csv",
             "external_survival_gse29609_summary.csv", "high_confidence_candidate_genes.csv",
             "external_survival_emtab1980.csv", "external_survival_gse29609.csv",
             "candidate_cox_bootstrap_summary.csv", "nested_cv_summary.csv",
             "candidate_clinical_composition_sensitivity.csv",
             "candidate_direct_tumor_purity_sensitivity.csv", "tumor_purity_coverage.csv",
             "hpa_candidate_cell_source_summary.csv", "nested_benchmark_repeat_metrics.csv")
metric <- function(x, key) as.numeric(x$value[match(key, x$metric)])
theme_review <- function() theme_minimal(base_size = 12, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, margin = margin(b = 12)),
        plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 12)),
        plot.margin = margin(15, 20, 15, 15))
save_pair <- function(p, stem, width, height) {
  ggsave(file.path(DIRS$figures, paste0(stem, ".png")), p, width = width,
         height = height, dpi = 300, bg = "white")
  ggsave(file.path(DIRS$figures, paste0(stem, ".pdf")), p, width = width,
         height = height, bg = "white")
}

# The six serial discovery filters stop at 23; the subsequent assessments branch.
funnel <- tab("evidence_funnel.csv")
stopifnot(identical(as.integer(funnel$count), c(8534L, 3323L, 1186L, 1117L, 538L, 23L)))
labels <- c("TCGA tumor vs normal differential expression",
            "GEO directions agree; nominal support in at least one",
            "Adjusted survival association (FDR < 0.05)",
            "Stage/grade direction and nominal significance",
            "Strict effect-size and consistency filters",
            "Stringent TCGA candidates (FDR < 0.01)")
funnel$y <- 8.2 - seq_len(nrow(funnel)) * 0.75
em_summary <- tab("external_survival_emtab1980_summary.csv")
gse_summary <- tab("external_survival_gse29609_summary.csv")
boot <- tab("candidate_cox_bootstrap_summary.csv")
cv <- tab("nested_cv_summary.csv")
purity <- tab("tumor_purity_coverage.csv")
branches <- data.frame(x = c(2.25, 6.75, 2.25, 6.75), y = c(2.5, 2.5, 1.1, 1.1),
  label = c(
    sprintf("External survival (unadjusted)\nE-MTAB: %d/%d, same-direction FDR\nGSE29609: %d/%d, same-direction nominal",
      metric(em_summary, "same_direction_fdr_candidates"), metric(em_summary, "platform_present_candidates"),
      metric(gse_summary, "same_direction_nominal_candidates"), metric(gse_summary, "platform_present_candidates")),
    sprintf("Fixed-candidate bootstrap\n%d/%d intervals exclude zero\nSelection was not repeated", sum(boot$ci_excludes_zero), nrow(boot)),
    sprintf("Held-out prediction: separate procedure\nNested selection, n = %d\nDelta C = %.4f [%.4f, %.4f]",
      cv$n_patients, cv$mean_delta_c, cv$patient_bootstrap_ci_low, cv$patient_bootstrap_ci_high),
    sprintf("Composition and cell-source assessment\nDirect purity: %d/%d retain direction + FDR\nMarker adjustment: %d/%d retain FDR\nNormal-tissue atlas; cell source remains uncertain",
      metric(purity, "purity_adjusted_fdr_candidates"), metric(purity, "revised_candidates_tested"),
      sum(tab("candidate_clinical_composition_sensitivity.csv")$composition_adjusted_fdr < 0.05), nrow(boot))))
p_funnel <- ggplot() +
  geom_rect(data = funnel, aes(xmin = 0.8, xmax = 8.2, ymin = y - 0.29, ymax = y + 0.29),
            fill = "#E9F1F3", color = "#A9C1C9", linewidth = 0.4) +
  geom_text(data = funnel, aes(x = 1, y = y, label = labels), hjust = 0, size = 3.65) +
  geom_text(data = funnel, aes(x = 8, y = y, label = format(count, big.mark = ",")),
            hjust = 1, size = 4.3, fontface = "bold", color = "#205D72") +
  annotate("text", x = 4.5, y = 3.25, label = "Downstream assessments are branches, not further candidate attrition",
           size = 3.6, fontface = "bold") +
  geom_rect(data = branches, aes(xmin = x - 2.1, xmax = x + 2.1, ymin = y - 0.58, ymax = y + 0.58),
            fill = "#F5F3EF", color = "#C9C2B7", linewidth = 0.4) +
  geom_text(data = branches, aes(x = x, y = y, label = label), size = 3.5, lineheight = 1.2) +
  coord_cartesian(xlim = c(0, 9), ylim = c(0.4, 8), clip = "off") +
  labs(title = "A stringent discovery funnel still leaves heterogeneous evidence",
       subtitle = "TCGA-KIRC discovery; GSE40435 and GSE53757 differential-expression replication",
       caption = "Boxes are schematic; counts are gene symbols. External denominators are platform-mapped candidates.\nBootstrap intervals condition on fixed selection. Prediction uses fold-specific selection; no clinical panel is validated.") +
  theme_review() + theme(axis.title = element_blank(), axis.text = element_blank(),
                         axis.ticks = element_blank(), panel.grid = element_blank())
save_pair(p_funnel, "review_evidence_funnel", 11, 9)

tcga <- tab("high_confidence_candidate_genes.csv")
em <- tab("external_survival_emtab1980.csv")
gse <- tab("external_survival_gse29609.csv")
effects <- bind_rows(
  tcga |> transmute(symbol, cohort = "TCGA-KIRC", estimate = main_log_hr,
                   low = log(main_hr_ci_low), high = log(main_hr_ci_high),
                   model = "Age + sex + stage + grade"),
  em |> filter(external_present) |> transmute(symbol, cohort = "E-MTAB-1980", estimate = external_log_hr,
                   low = log(external_hr_ci_low), high = log(external_hr_ci_high), model = "Expression only"),
  gse |> filter(external_present) |> transmute(symbol, cohort = "GSE29609", estimate = external_log_hr,
                   low = log(external_hr_ci_low), high = log(external_hr_ci_high), model = "Expression only")) |>
  mutate(row = paste0(cohort, "\n", model),
         row = factor(row, levels = rev(c("TCGA-KIRC\nAge + sex + stage + grade",
                    "E-MTAB-1980\nExpression only", "GSE29609\nExpression only"))))
forest_plot <- function(symbol, title, subtitle) {
  d <- effects |> filter(.data$symbol == .env$symbol)
  stopifnot(nrow(d) == 3, all(is.finite(d$low)), all(d$low < d$high))
  ggplot(d, aes(x = exp(estimate), y = row, color = cohort)) +
    geom_vline(xintercept = 1, color = "#7A7A7A", linetype = "dashed") +
    geom_errorbar(aes(xmin = exp(low), xmax = exp(high)), orientation = "y", width = 0.12, linewidth = 0.7) +
    geom_point(size = 3) + scale_x_log10() +
    scale_color_manual(values = c("TCGA-KIRC" = "#205D72", "E-MTAB-1980" = "#3E6B48", "GSE29609" = "#B85C38")) +
    labs(title = title, subtitle = subtitle, x = "Hazard ratio per cohort-specific SD of expression (95% CI)", y = NULL,
      caption = "TCGA: n = 517, 170 deaths; E-MTAB: n = 101, 23 deaths; GSE29609: n = 39, 17 deaths.\nExternal plotted models are unadjusted; adjusted sensitivity is reported in the evidence ledger.") +
    theme_review() + theme(legend.position = "none")
}
save_pair(forest_plot("ACADM", "ACADM: replication in E-MTAB does not establish uniform transportability",
  "Lower hazard in TCGA and E-MTAB; the smaller GSE29609 cohort is inconclusive and points oppositely"),
  "review_replicated_association", 10.5, 4.8)

# A two-panel figure preserves the distinction between cohort contradiction and
# marker-score attenuation. The stored marker-score table has no CI columns.
comp <- tab("candidate_clinical_composition_sensitivity.csv") |> filter(symbol == "DDC")
direct <- tab("candidate_direct_tumor_purity_sensitivity.csv") |> filter(symbol == "DDC")
stopifnot(nrow(comp) == 1, nrow(direct) == 1)
ddc <- effects |> filter(symbol == "DDC") |>
  transmute(panel = "A  Cross-cohort estimates with 95% CIs", row = as.character(row),
            estimate, low, high)
composition <- data.frame(panel = "B  TCGA adjustment sensitivity (coefficients only)",
  row = c("Clinical model\nn = 517", "Clinical + marker scores\nn = 517",
          "Matched clinical model\nn = 516", "Clinical + direct purity\nn = 516"),
  estimate = c(comp$gene_log_hr, comp$composition_adjusted_log_hr,
               direct$matched_baseline_log_hr, direct$gene_log_hr), low = NA_real_, high = NA_real_)
d <- bind_rows(ddc, composition)
d$row <- factor(d$row, levels = rev(c(as.character(effects$row[effects$symbol == "DDC"]), composition$row)))
p_ddc <- ggplot(d, aes(estimate, row)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7A7A7A") +
  geom_errorbar(data = d |> filter(is.finite(low)), aes(xmin = low, xmax = high),
                orientation = "y", width = 0.1, linewidth = 0.7, color = "#B85C38") +
  geom_point(size = 3, color = "#B85C38") +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  labs(title = "DDC: direction reverses externally and marker adjustment attenuates the association",
       subtitle = sprintf("Marker-score adjusted log HR = %.3f (FDR %.3f); direct purity changes the matched coefficient by %.3f",
         comp$composition_adjusted_log_hr, comp$composition_adjusted_fdr,
         abs(direct$gene_log_hr - direct$matched_baseline_log_hr)),
       x = "Log hazard ratio per cohort-specific SD of expression", y = NULL,
       caption = "A: TCGA adjusted model; external expression-only models. GSE29609: n = 39, 17 deaths.\nB: Marker-score coefficients have no stored CIs; dots do not imply equal precision. Direct purity and marker scores address different threats.\nNormal-tissue atlas cell-source context is descriptive and cannot establish the source of tumor expression.") +
  theme_review() + theme(strip.text = element_text(hjust = 0, face = "bold"),
                         panel.spacing = grid::unit(1, "lines"))
save_pair(p_ddc, "review_contradictory_case", 11.5, 8.2)

# Refresh the maintained benchmark filename from recorded fold-wise fits.
# This display-only step does not change checkpoint signatures or refit models.
benchmark <- tab("nested_benchmark_repeat_metrics.csv") |>
  filter(strategy != "clinical") |>
  mutate(strategy = factor(strategy,
    levels = c("hydra", "survival_only", "de_only", "ridge_eligible", "matched_control"),
    labels = c("HYDRA", "Survival only", "DE only", "Ridge on eligible DEGs", "Matched control")))
stopifnot(!anyNA(benchmark$strategy), n_distinct(benchmark$strategy) == 5L)
p_benchmark <- ggplot(benchmark, aes(strategy, delta_c)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, color = "#9b3d2a", linetype = "dashed", linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
             color = "#216b75", size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 10, color = "#173f46") +
  labs(x = NULL, y = "Held-out concordance minus clinical",
       title = "Nested selection benchmark",
       subtitle = "Same outer splits; dashed line is the project 0.01 increment. Points are repeats.",
       caption = "Display regenerated from nested_benchmark_repeat_metrics.csv; model fitting remains in stage 31.\nRepeat-level variation is not an independent-cohort confidence interval.") +
  theme_review() + theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggsave(file.path(DIRS$figures, "nested_selection_benchmark.png"), p_benchmark,
       width = 8.4, height = 4.8, dpi = 180, bg = "white")

# A small provenance sidecar makes the figure dependencies and runtime explicit.
write_csv(data.frame(source = file.path(DIRS$tables, sources),
                     md5 = unname(tools::md5sum(file.path(DIRS$tables, sources)))),
          file.path(DIRS$tables, "review_figure_sources.csv"))
writeLines(c("Command: Rscript analysis/40_review_figures.R",
             "All estimates are stored model outputs; benchmark point jitter uses seed 1.",
             capture.output(sessionInfo())), "environment/review_figures_sessionInfo.txt")
message("Generated three review figures in PNG/PDF and refreshed the maintained benchmark PNG, with provenance.")
