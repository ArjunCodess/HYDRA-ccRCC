source("analysis/00_config.R")
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
})

cv <- read_csv(file.path(DIRS$tables, "nested_cv_repeat_metrics.csv"), show_col_types = FALSE)
summary <- read_csv(file.path(DIRS$tables, "nested_cv_summary.csv"), show_col_types = FALSE)
funnel <- read_csv(file.path(DIRS$tables, "funnel_external_summary.csv"), show_col_types = FALSE)
rule_labels <- c(
  complete_rule = "Complete",
  survival_only = "Survival only",
  de_only = "DE only",
  leave_out_gse40435 = "Leave GSE40435 out",
  leave_out_gse53757 = "Leave GSE53757 out",
  expression_matched_control = "Expression-matched control"
)
stopifnot(setequal(funnel$rule, names(rule_labels)))
funnel$rule <- factor(funnel$rule, levels = names(rule_labels),
                      labels = unname(rule_labels))

p1 <- ggplot(cv, aes(repeat_id, delta_c)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, color = "#9b3d2a", linetype = "dashed", linewidth = 0.5) +
  geom_point(color = "#216b75", size = 2.2) +
  geom_hline(yintercept = summary$mean_delta_c, color = "#216b75", linewidth = 0.7) +
  scale_x_continuous(breaks = seq_len(nrow(cv))) +
  labs(x = "Outer-CV repeat", y = "Selected-gene minus clinical concordance",
       title = "Selection-aware held-out increment",
       subtitle = "Dashed line: prespecified 0.01 increment threshold") +
  theme_minimal(base_size = 12)
ggsave(file.path(DIRS$figures, "nested_cv_increment.png"), p1,
       width = 8, height = 4.5, dpi = 180)

p2 <- ggplot(funnel, aes(rule, directional_rate, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("E-MTAB-1980" = "#216b75", "GSE29609" = "#b77443")) +
  labs(x = NULL, y = "Same TCGA hazard direction among mapped genes",
       title = "Exploratory external comparison at matched list size",
       fill = "Cohort") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 23, hjust = 1))
ggsave(file.path(DIRS$figures, "funnel_external_comparison.png"), p2,
       width = 9, height = 4.7, dpi = 180)
message("Audit figures regenerated.")
