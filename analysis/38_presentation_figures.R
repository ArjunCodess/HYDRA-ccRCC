# Presentation figures from saved tables. No new model fits.

source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

tab <- function(name) read_csv(file.path(DIRS$tables, name), show_col_types = FALSE)
value_of <- function(data, metric) {
  hit <- data$value[data$metric == metric]
  if (length(hit) != 1L || !is.finite(hit)) stop("Missing metric: ", metric)
  hit
}
expect_near <- function(x, target, tol = 5e-4, label = "value") {
  if (!is.finite(x) || abs(x - target) > tol) {
    stop(label, " is ", x, "; expected ", target)
  }
  invisible(x)
}
signed4 <- function(x) sprintf("%+.4f", x)
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

ink <- "#1c1c1c"
muted <- "#5c5c5c"
rule_col <- "#d9d9d9"
hydra_col <- "#1B4F72"
ridge_col <- "#8C2F39"
neutral <- "#8d8d8d"
paper <- "#f7f4ef"
chip <- "#efeae3"

arithmetic <- tab("tcga_aliquot_arithmetic.csv")
funnel <- tab("candidate_summary.csv")
repro <- tab("reproducibility_summary.csv")
bench <- tab("nested_benchmark_summary.csv")
clearcode <- tab("published_signature_summary.csv") |>
  filter(cohort == "TCGA-KIRC", arm == "clinical_plus_clearcode34")
gse_s <- tab("external_survival_gse29609_summary.csv")
em_s <- tab("external_survival_emtab1980_summary.csv")
unpaired <- tab("gse53757_unpaired_summary.csv")
aliquot <- tab("aliquot_sensitivity_summary.csv")
dispersion <- tab("tracerx_direction_dispersion.csv")
evidence <- tab("candidate_evidence_matrix.csv")

n_aliquot <- value_of(arithmetic, "tumor_aliquots")
n_patient <- value_of(arithmetic, "tumor_patients")
n_pair <- value_of(arithmetic, "normal_patients")
n_mapped <- value_of(repro, "tcga_significant")
n_repro <- value_of(funnel, "reproducible_deg")
n_surv <- value_of(funnel, "main_stage_grade_complete_prognostic")
n_strict <- value_of(funnel, "strict_candidate")
n_high <- value_of(funnel, "high_confidence_candidate")
n_nested <- 517L

expect_near(n_aliquot, 541, 0, "aliquots")
expect_near(n_patient, 533, 0, "patients")
expect_near(n_pair, 72, 0, "pairs")
expect_near(n_mapped, 8534, 0, "mapped DEGs")
expect_near(n_repro, 3323, 0, "reproducible")
expect_near(n_surv, 1186, 0, "survival")
expect_near(n_strict, 538, 0, "strict")
expect_near(n_high, 23, 0, "high-confidence")
if (nrow(evidence) != 23L) stop("Evidence matrix does not have 23 genes.")

arm_of <- function(strategy) {
  row <- bench[bench$strategy == strategy, ]
  if (nrow(row) != 1L) stop("Missing benchmark arm: ", strategy)
  row
}
hydra <- arm_of("hydra")
survival_only <- arm_of("survival_only")
de_only <- arm_of("de_only")
matched <- arm_of("matched_control")
ridge <- arm_of("ridge_eligible")
expect_near(hydra$mean_delta_c, 0.0039, label = "HYDRA")
expect_near(survival_only$mean_delta_c, -0.0011, label = "survival-only")
expect_near(de_only$mean_delta_c, -0.0015, label = "DE-only")
expect_near(matched$mean_delta_c, 0.0019, label = "matched")
expect_near(clearcode$delta_c, 0.0129, label = "ClearCode34")
expect_near(ridge$mean_delta_c, 0.0257, label = "ridge")
expect_near(value_of(gse_s, "same_direction_candidates"), 5, 0, "GSE29609 same")
expect_near(value_of(gse_s, "platform_present_candidates"), 21, 0, "GSE29609 mapped")
expect_near(value_of(em_s, "same_direction_candidates"), 21, 0, "E-MTAB same")
expect_near(value_of(em_s, "platform_present_candidates"), 22, 0, "E-MTAB mapped")
expect_near(value_of(em_s, "strict_external_support_candidates"), 12, 0, "E-MTAB strict")
expect_near(value_of(unpaired, "high_confidence_still_reproducible"), 23, 0, "unpaired retained")
expect_near(value_of(unpaired, "reproducible_jaccard"), 0.998, 0.0005, "Jaccard")

arms <- tibble(
  arm = c("HYDRA", "Survival-only", "DE-only", "Matched control", "ClearCode34", "Ridge"),
  delta = c(hydra$mean_delta_c, survival_only$mean_delta_c, de_only$mean_delta_c,
            matched$mean_delta_c, clearcode$delta_c, ridge$mean_delta_c),
  lo = c(hydra$patient_bootstrap_ci_low, survival_only$patient_bootstrap_ci_low,
         de_only$patient_bootstrap_ci_low, matched$patient_bootstrap_ci_low,
         clearcode$delta_c_ci_low, ridge$patient_bootstrap_ci_low),
  hi = c(hydra$patient_bootstrap_ci_high, survival_only$patient_bootstrap_ci_high,
         de_only$patient_bootstrap_ci_high, matched$patient_bootstrap_ci_high,
         clearcode$delta_c_ci_high, ridge$patient_bootstrap_ci_high),
  emphasis = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

draw_master <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  steps <- tibble(
    n = c(n_aliquot, n_patient, n_pair, n_mapped, n_repro, n_surv, n_strict, n_high, n_nested),
    label = c("aliquots", "patients", "pairs", "mapped", "reproducible", "survival", "strict", "candidates", "nested")
  )
  n_step <- nrow(steps)
  xs <- seq(0.055, 0.945, length.out = n_step)
  grid::grid.text("541  to  the held-out benchmark", x = 0.04, y = 0.955,
                  just = c("left", "center"),
                  gp = grid::gpar(cex = 1.15, fontface = "bold", col = ink))
  for (i in seq_len(n_step)) {
    if (i < n_step) {
      grid::grid.lines(
        x = grid::unit(c(xs[i] + 0.028, xs[i + 1] - 0.028), "npc"),
        y = grid::unit(c(0.78, 0.78), "npc"),
        arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed"),
        gp = grid::gpar(col = "#b9b3aa", fill = "#b9b3aa", lwd = 1.4)
      )
    }
    fill <- if (steps$label[i] == "nested") "#E7EEF2" else chip
    border <- if (steps$label[i] == "nested") hydra_col else "#d5cfc6"
    grid::grid.roundrect(
      x = xs[i], y = 0.78, width = 0.09, height = 0.22,
      r = grid::unit(0.012, "snpc"),
      gp = grid::gpar(fill = fill, col = border, lwd = 1.3)
    )
    grid::grid.text(comma(steps$n[i]), x = xs[i], y = 0.80,
                    gp = grid::gpar(cex = 1.35, fontface = "bold", col = ink))
    grid::grid.text(steps$label[i], x = xs[i], y = 0.70,
                    gp = grid::gpar(cex = 0.72, col = muted))
  }

  plot_vp <- grid::viewport(x = 0.50, y = 0.36, width = 0.94, height = 0.52)
  grid::pushViewport(plot_vp)
  xlim <- c(-0.020, 0.050)
  ylim <- c(0.45, nrow(arms) + 0.55)
  x_at <- function(v) (v - xlim[1]) / diff(xlim)
  y_at <- function(v) (v - ylim[1]) / diff(ylim)
  grid::grid.rect(gp = grid::gpar(fill = paper, col = NA))
  grid::grid.lines(x = x_at(c(0, 0)), y = c(0.08, 0.92), gp = grid::gpar(col = ink, lwd = 0.7))
  grid::grid.lines(x = x_at(c(0.01, 0.01)), y = c(0.08, 0.92),
                   gp = grid::gpar(col = ridge_col, lwd = 0.8, lty = "dashed"))
  grid::grid.text("+0.01", x = x_at(0.01), y = 0.955, gp = grid::gpar(cex = 0.68, col = ridge_col))
  for (i in seq_len(nrow(arms))) {
    row <- arms[i, ]
    yy <- y_at(nrow(arms) - i + 1)
    col <- if (row$arm == "HYDRA") hydra_col else if (row$arm == "Ridge") ridge_col else neutral
    lwd <- if (row$emphasis) 2.4 else 1.3
    grid::grid.lines(
      x = x_at(c(row$lo, row$hi)), y = c(yy, yy),
      gp = grid::gpar(col = col, lwd = lwd, lineend = "butt")
    )
    grid::grid.points(x = x_at(row$delta), y = yy, pch = 16, size = grid::unit(0.55, "char"),
                      gp = grid::gpar(col = col))
    grid::grid.text(row$arm, x = 0.012, y = yy, just = "left",
                    gp = grid::gpar(cex = if (row$emphasis) 1.05 else 0.9,
                                    fontface = if (row$emphasis) "bold" else "plain",
                                    col = if (row$emphasis) col else muted))
    grid::grid.text(signed4(row$delta), x = 0.985, y = yy, just = "right",
                    gp = grid::gpar(cex = if (row$emphasis) 1.25 else 0.9,
                                    fontface = if (row$emphasis) "bold" else "plain",
                                    col = col))
  }
  grid::popViewport()
  grid::grid.text("Same 517 patients and splits. Change in concordance versus the clinical model.",
                  x = 0.04, y = 0.035, just = "left",
                  gp = grid::gpar(cex = 0.78, col = muted))
}

save_master <- function() {
  dir.create("paper/figures", recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file.path(DIRS$figures, "master_funnel_benchmark.png"),
                 width = 16.4, height = 7.1, units = "in", res = 220)
  draw_master()
  grDevices::dev.off()
  grDevices::pdf(file.path("paper/figures", "hydra_evidence_overview.pdf"),
                 width = 16.4, height = 7.1)
  draw_master()
  grDevices::dev.off()
  grDevices::svg(file.path("paper/figures", "hydra_evidence_overview.svg"),
                 width = 16.4, height = 7.1)
  draw_master()
  grDevices::dev.off()
}

draw_external <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("External direction", x = 0.05, y = 0.90, just = "left",
                  gp = grid::gpar(cex = 1.3, fontface = "bold", col = ink))
  blocks <- list(
    list(x = 0.27, num = "5/21", sub = "GSE29609\nsame direction", fill = "#f3e6df", col = "#8C2F39"),
    list(x = 0.73, num = "21/22", sub = "E-MTAB-1980\nsame direction", fill = "#e5efe8", col = "#2E5A3C")
  )
  for (block in blocks) {
    grid::grid.roundrect(x = block$x, y = 0.58, width = 0.38, height = 0.42,
                         r = grid::unit(0.02, "snpc"),
                         gp = grid::gpar(fill = block$fill, col = block$col, lwd = 1.6))
    grid::grid.text(block$num, x = block$x, y = 0.62,
                    gp = grid::gpar(cex = 3.4, fontface = "bold", col = ink))
    grid::grid.text(block$sub, x = block$x, y = 0.45,
                    gp = grid::gpar(cex = 0.95, col = muted))
  }
  grid::grid.text("12 strict", x = 0.73, y = 0.72,
                  gp = grid::gpar(cex = 1.05, fontface = "bold", col = "#2E5A3C"))
  grid::grid.roundrect(x = 0.50, y = 0.16, width = 0.72, height = 0.16,
                       r = grid::unit(0.015, "snpc"),
                       gp = grid::gpar(fill = "#F6E4DC", col = "#8C2F39", lwd = 1.2))
  grid::grid.text("DDC  and  TCIRG1   reverse in GSE29609", x = 0.50, y = 0.16,
                  gp = grid::gpar(cex = 1.15, fontface = "bold", col = "#8C2F39"))
}

draw_cohort <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("533 patients  to  517 nested", x = 0.04, y = 0.92, just = "left",
                  gp = grid::gpar(cex = 1.2, fontface = "bold", col = ink))
  boxes <- list(
    list(x = 0.16, y = 0.62, w = 0.24, h = 0.28, title = "533", sub = "patients\n175 deaths", fill = chip, col = "#b9b3aa"),
    list(x = 0.84, y = 0.62, w = 0.24, h = 0.28, title = "517", sub = "complete cases\n170 deaths", fill = "#E7EEF2", col = hydra_col)
  )
  for (box in boxes) {
    grid::grid.roundrect(x = box$x, y = box$y, width = box$w, height = box$h,
                         r = grid::unit(0.02, "snpc"),
                         gp = grid::gpar(fill = box$fill, col = box$col, lwd = 1.5))
    grid::grid.text(box$title, x = box$x, y = box$y + 0.04,
                    gp = grid::gpar(cex = 2.2, fontface = "bold", col = ink))
    grid::grid.text(box$sub, x = box$x, y = box$y - 0.06,
                    gp = grid::gpar(cex = 0.85, col = muted))
  }
  drops <- tibble(
    x = c(0.34, 0.46, 0.58, 0.70),
    n = c("8", "3", "1", "4"),
    why = c("no grade", "no stage", "no age", "time = 0")
  )
  for (i in seq_len(nrow(drops))) {
    grid::grid.roundrect(x = drops$x[i], y = 0.28, width = 0.11, height = 0.22,
                         r = grid::unit(0.012, "snpc"),
                         gp = grid::gpar(fill = "#F6E4DC", col = "#8C2F39", lwd = 1))
    grid::grid.text(drops$n[i], x = drops$x[i], y = 0.31,
                    gp = grid::gpar(cex = 1.35, fontface = "bold", col = "#8C2F39"))
    grid::grid.text(drops$why[i], x = drops$x[i], y = 0.22,
                    gp = grid::gpar(cex = 0.68, col = muted))
  }
  grid::grid.text("16 excluded, one field each.  5 of the 16 died.",
                  x = 0.50, y = 0.08, gp = grid::gpar(cex = 0.9, col = ink))
}

draw_pairing <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("GSE53757", x = 0.05, y = 0.88, just = "left",
                  gp = grid::gpar(cex = 1.2, fontface = "bold", col = ink))
  panels <- list(
    list(x = 0.27, num = "23/23", sub = "high-confidence genes\nretained unpaired"),
    list(x = 0.73, num = "0.998", sub = "Jaccard overlap of\nthe reproducible lists")
  )
  for (panel in panels) {
    grid::grid.roundrect(x = panel$x, y = 0.48, width = 0.40, height = 0.52,
                         r = grid::unit(0.02, "snpc"),
                         gp = grid::gpar(fill = chip, col = "#d5cfc6", lwd = 1.3))
    grid::grid.text(panel$num, x = panel$x, y = 0.54,
                    gp = grid::gpar(cex = 2.8, fontface = "bold", col = ink))
    grid::grid.text(panel$sub, x = panel$x, y = 0.34,
                    gp = grid::gpar(cex = 0.9, col = muted))
  }
}

draw_aliquot <- function() {
  primary <- 23L
  rules <- tibble(
    rule = c("deepest", "drop_multi", "first_barcode", "sum_counts"),
    label = c("Deepest library", "Drop duplicated\npatients", "First barcode", "Summed counts"),
    n = c(primary,
          aliquot$n_high_confidence[aliquot$rule == "drop_multi"],
          aliquot$n_high_confidence[aliquot$rule == "first_barcode"],
          aliquot$n_high_confidence[aliquot$rule == "sum_counts"])
  )
  if (any(rules$n != c(23, 19, 20, 19))) stop("Aliquot high-confidence counts changed.")
  lost <- strsplit(aliquot$lost, ";")
  lost_all <- Reduce(intersect, lost)
  if (!setequal(lost_all, c("GRAMD1A", "IFFO1", "LTB4R", "RBM47"))) {
    stop("Genes lost under every aliquot rule changed.")
  }
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("Aliquot rule", x = 0.05, y = 0.90, just = "left",
                  gp = grid::gpar(cex = 1.15, fontface = "bold", col = ink))
  xs <- seq(0.16, 0.84, length.out = 4)
  for (i in seq_len(4)) {
    primary_rule <- i == 1L
    grid::grid.roundrect(x = xs[i], y = 0.58, width = 0.18, height = 0.40,
                         r = grid::unit(0.015, "snpc"),
                         gp = grid::gpar(fill = if (primary_rule) "#E7EEF2" else chip,
                                         col = if (primary_rule) hydra_col else "#d5cfc6", lwd = 1.4))
    grid::grid.text(rules$n[i], x = xs[i], y = 0.64,
                    gp = grid::gpar(cex = 2.1, fontface = "bold", col = ink))
    grid::grid.text(rules$label[i], x = xs[i], y = 0.48,
                    gp = grid::gpar(cex = 0.72, col = muted))
  }
  grid::grid.text("Lost under every alternate rule:  GRAMD1A, IFFO1, LTB4R, RBM47",
                  x = 0.50, y = 0.16, gp = grid::gpar(cex = 0.95, fontface = "bold", col = "#8C2F39"))
}

draw_tracerx <- function() {
  repeats <- tab("tracerx_one_region_cox_repeats.csv") |>
    filter(scenario %in% c("fixed_subset_regions", "size_matched_39")) |>
    group_by(scenario, repeat_id) |>
    summarise(agreement = mean(same_tcga_direction), .groups = "drop")
  fixed <- dispersion$median[dispersion$scenario == "fixed_subset_regions"]
  changing <- dispersion$median[dispersion$scenario == "size_matched_39"]
  fixed_lo <- dispersion$iqr_low[dispersion$scenario == "fixed_subset_regions"]
  fixed_hi <- dispersion$iqr_high[dispersion$scenario == "fixed_subset_regions"]
  ch_lo <- dispersion$iqr_low[dispersion$scenario == "size_matched_39"]
  ch_hi <- dispersion$iqr_high[dispersion$scenario == "size_matched_39"]
  expect_near(fixed, 0.985, 0.0005, "TRACERx fixed median")
  expect_near(changing, 0.863, 0.0005, "TRACERx changing median")
  labels <- c(
    fixed_subset_regions = "98.5%   [76.6–100]",
    size_matched_39 = "86.3%   [77.8–93.1]"
  )
  subtitles <- c(
    fixed_subset_regions = "39 patients fixed. Region changes. Across genes.",
    size_matched_39 = "39 patients redrawn. Region changes. Across genes."
  )
  repeats$panel <- labels[repeats$scenario]
  repeats$panel <- factor(repeats$panel, levels = unname(labels))
  ggplot2::ggplot(repeats, ggplot2::aes(agreement)) +
    ggplot2::geom_histogram(bins = 28, fill = hydra_col, color = "white", linewidth = 0.15) +
    ggplot2::facet_wrap(~ panel, ncol = 1) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Candidates keeping the TCGA hazard direction, per draw", y = "Draws",
                  title = "TRACERx direction",
                  subtitle = "Panel labels are the median and interquartile range across genes.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 13),
      panel.grid.minor = ggplot2::element_blank()
    )
}

fmt_num <- function(x, digits = 2) {
  if (!is.finite(x)) return("---")
  formatC(x, format = "f", digits = digits)
}
direction_label <- function(present, same) {
  ifelse(!present, "absent", ifelse(same, "same", "opposite"))
}
evidence_tex <- function() {
  lost_all <- c("GRAMD1A", "IFFO1", "LTB4R", "RBM47")
  rows <- evidence |>
    arrange(desc(evidence_score), symbol) |>
    mutate(
      gene = ifelse(symbol %in% c("DDC", "TCIRG1"), paste0("\\textbf{", symbol, "}"), symbol),
      gene = ifelse(symbol %in% lost_all, paste0(gene, "$^{*}$"), gene),
      gse_dir = direction_label(gse29609_present, gse29609_same_direction),
      em_dir = direction_label(emtab_present, emtab_same_direction)
    )
  if (sum(rows$gse_dir == "same") != 5L || sum(rows$em_dir == "same") != 21L) {
    stop("Evidence-matrix directions do not match 5/21 and 21/22.")
  }
  if (!all(rows$gse_dir[rows$symbol %in% c("DDC", "TCIRG1")] == "opposite")) {
    stop("DDC and TCIRG1 are not marked as GSE29609 reversals.")
  }
  body <- apply(rows, 1, function(row) {
    paste(
      row[["gene"]],
      fmt_num(as.numeric(row[["tcga_log2fc"]])),
      fmt_num(as.numeric(row[["main_hr"]])),
      ifelse(isTRUE(as.logical(row[["gse40435_nominal"]])), "yes", "no"),
      ifelse(isTRUE(as.logical(row[["gse53757_nominal"]])), "yes", "no"),
      row[["gse_dir"]],
      row[["em_dir"]],
      ifelse(isTRUE(as.logical(row[["marker_fdr_support"]])), "yes", "no"),
      ifelse(isTRUE(as.logical(row[["purity_fdr_support"]])), "yes", "no"),
      row[["nested_folds"]],
      sprintf("%.0f\\%%", 100 * as.numeric(row[["tracerx_fixed_direction"]])),
      ifelse(isTRUE(as.logical(row[["retained_under_every_aliquot_rule"]])), "yes", "no"),
      sep = " & "
    )
  })
  header <- paste(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\caption{Evidence matrix for the 23 high-confidence candidates. Every row uses the same columns. Bold marks the GSE29609 reversals DDC and TCIRG1. An asterisk marks a gene lost under every alternate aliquot rule. Direction is relative to the TCGA hazard. Nested folds count how often the gene was the selected gene. TRACERx is the fraction of fixed-patient region draws that kept the TCGA direction.}",
    "\\label{tab:evidence}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrcccccccrr}",
    "\\toprule",
    "Gene & TCGA log$_2$FC & TCGA HR & GSE40435 & GSE53757 & GSE29609 & E-MTAB-1980 & Marker & Purity & Folds & TRACERx & Aliquot \\\\",
    "\\midrule",
    paste0(body, " \\\\", collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}",
    sep = "\n"
  )
  writeLines(header, "paper/evidence_matrix.tex")
}

save_png <- function(path, width, height, draw) {
  grDevices::png(path, width = width, height = height, units = "in", res = 220)
  draw()
  grDevices::dev.off()
}

save_master()
save_png(file.path(DIRS$figures, "external_direction_comparison.png"), 8.6, 4.6, draw_external)
save_png(file.path(DIRS$figures, "cohort_flow.png"), 9.2, 4.2, draw_cohort)
save_png(file.path(DIRS$figures, "gse53757_pairing_sensitivity.png"), 8.4, 4.2, draw_pairing)
save_png(file.path(DIRS$figures, "aliquot_sensitivity.png"), 9.0, 4.2, draw_aliquot)
ggplot2::ggsave(file.path(DIRS$figures, "tracerx_sampling_distributions.png"),
                draw_tracerx(), width = 7.6, height = 5.4, dpi = 220)
evidence_tex()
message("Presentation figures written.")
