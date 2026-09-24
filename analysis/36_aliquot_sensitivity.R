# Aliquot-rule sensitivity. Does not replace the primary deepest-library shortlist.
# One rule: HYDRA_ALIQUOT_RULE=drop_multi|first_barcode|sum_counts
# All three when the variable is unset.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(tibble)
  library(survival)
})

se <- read_required_rds(FILES$tcga_se)
raw_counts <- assay(se, "unstranded")
raw_coldata <- as.data.frame(colData(se)) |>
  tibble::rownames_to_column("sample_barcode")
raw_coldata$patient_barcode <- substr(raw_coldata$sample_barcode, 1, 12)
raw_coldata$library_depth <- as.numeric(colSums(raw_counts[, raw_coldata$sample_barcode, drop = FALSE]))

tumors <- raw_coldata |> filter(sample_type == "Primary Tumor")
multi <- tumors |> count(patient_barcode) |> filter(n > 1L) |> pull(patient_barcode)
audit <- tumors |>
  filter(patient_barcode %in% multi) |>
  arrange(patient_barcode, desc(library_depth), sample_barcode) |>
  group_by(patient_barcode) |>
  mutate(selected_by_deepest = row_number() == 1L) |>
  ungroup() |>
  select(patient_barcode, sample_barcode, sample_type, library_depth, selected_by_deepest)
write_csv_atomic(audit, file.path(DIRS$tables, "tcga_aliquot_selection_audit.csv"))

normals <- raw_coldata |> filter(sample_type == "Solid Tissue Normal")
write_csv_atomic(tibble(
  metric = c("tumor_aliquots", "tumor_patients", "excess_tumor_aliquots",
             "patients_with_three_tumor_aliquots", "normal_aliquots",
             "normal_patients", "patients_with_multiple_normals"),
  value = c(nrow(tumors), n_distinct(tumors$patient_barcode),
            nrow(tumors) - n_distinct(tumors$patient_barcode),
            sum(table(tumors$patient_barcode) == 3L),
            nrow(normals), n_distinct(normals$patient_barcode),
            sum(table(normals$patient_barcode) > 1L))
), file.path(DIRS$tables, "tcga_aliquot_arithmetic.csv"))

rules <- Sys.getenv("HYDRA_ALIQUOT_RULE", "")
rules <- if (nzchar(rules)) strsplit(rules, ",", fixed = TRUE)[[1]] else {
  c("drop_multi", "first_barcode", "sum_counts")
}

clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE) |>
  transmute(patient_barcode = submitter_id, os_time, os_event,
            age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
            sex = gender,
            stage = normalize_stage(ajcc_pathologic_stage),
            grade_raw = normalize_grade(tumor_grade),
            grade = case_when(grade_raw %in% c("G1", "G2") ~ "Low grade",
                              grade_raw %in% c("G3", "G4") ~ "High grade",
                              TRUE ~ NA_character_))
primary <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
                    show_col_types = FALSE) |> pull(symbol)
geo_a <- read_csv(file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"), show_col_types = FALSE)
geo_b <- read_csv(file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"), show_col_types = FALSE)
out_dir <- file.path(DIRS$tables, "aliquot_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_rule <- function(rule) {
  message("Aliquot rule: ", rule)
  picked <- apply_aliquot_rule(raw_coldata, raw_counts, rule)
  coldata <- as.data.frame(picked$coldata)
  counts <- picked$counts
  coldata$condition <- factor(ifelse(coldata$sample_type == "Primary Tumor", "TP", "NT"),
                              levels = c("NT", "TP"))
  rownames(coldata) <- coldata$sample_barcode
  counts <- counts[, coldata$sample_barcode, drop = FALSE]
  keep <- rowSums(counts >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
  dds <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], coldata, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "TP", "NT")))
  res$gene_id <- rownames(res)
  vsd <- assay(vst(dds, blind = TRUE))
  symbols <- AnnotationDbi::select(org.Hs.eg.db,
    keys = unique(sub("\\..*$", "", res$gene_id)),
    keytype = "ENSEMBL", columns = "SYMBOL") |>
    filter(!is.na(SYMBOL)) |> arrange(ENSEMBL, SYMBOL) |> distinct(ENSEMBL, .keep_all = TRUE)
  deg <- res |>
    mutate(ensembl = sub("\\..*$", "", gene_id),
           tcga_significant = !is.na(padj) & padj < THRESHOLDS$deg_fdr &
             abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc,
           tcga_direction = sign(log2FoldChange)) |>
    left_join(symbols, by = c("ensembl" = "ENSEMBL")) |>
    filter(!is.na(SYMBOL)) |>
    group_by(SYMBOL) |> arrange(padj, .by_group = TRUE) |> slice_head(n = 1) |> ungroup() |>
    transmute(symbol = SYMBOL, gene_id, tcga_log2fc = log2FoldChange, tcga_padj = padj,
              tcga_significant, tcga_direction)
  joined <- deg |>
    inner_join(geo_a |> transmute(symbol, gse40435_log2fc = log2FoldChange, gse40435_pvalue = pvalue), by = "symbol") |>
    inner_join(geo_b |> transmute(symbol, gse53757_log2fc = log2FoldChange, gse53757_pvalue = pvalue), by = "symbol") |>
    mutate(reproducible_deg = tcga_significant &
             sign(gse40435_log2fc) == tcga_direction & sign(gse53757_log2fc) == tcga_direction &
             (gse40435_pvalue < 0.05 | gse53757_pvalue < 0.05)) |>
    filter(reproducible_deg, gene_id %in% rownames(vsd))

  tumors_only <- coldata |>
    filter(condition == "TP") |>
    inner_join(clinical, by = "patient_barcode") |>
    filter(is.finite(os_time), os_time > 0, os_event %in% c(0L, 1L),
           is.finite(age), !is.na(sex), !is.na(stage), !is.na(grade))
  expr_mat <- vsd[, tumors_only$sample_barcode, drop = FALSE]
  specs <- list(stage_grade_complete = c("age", "sex", "stage", "grade"),
                stage_complete = c("age", "sex", "stage"),
                grade_complete = c("age", "sex", "grade"))
  one_model <- function(gene_id, covars) {
    dat <- tumors_only
    dat$expr <- as.numeric(scale(expr_mat[gene_id, ]))
    dat$sex <- factor(dat$sex)
    dat$stage <- factor(dat$stage)
    dat$grade <- factor(dat$grade)
    fit <- tryCatch(suppressWarnings(coxph(
      as.formula(paste("Surv(os_time, os_event) ~ expr +", paste(covars, collapse = " + "))),
      data = dat)), error = function(e) NULL)
    if (is.null(fit) || !"expr" %in% names(coef(fit))) return(c(beta = NA_real_, p = NA_real_))
    c(beta = unname(coef(fit)["expr"]), p = summary(fit)$coefficients["expr", "Pr(>|z|)"])
  }
  main <- t(vapply(joined$gene_id, one_model, numeric(2), covars = specs$stage_grade_complete))
  joined$main_beta <- main[, 1]
  joined$main_p <- main[, 2]
  joined$main_fdr <- p.adjust(joined$main_p, method = "BH")
  stage <- t(vapply(joined$gene_id, one_model, numeric(2), covars = specs$stage_complete))
  grade <- t(vapply(joined$gene_id, one_model, numeric(2), covars = specs$grade_complete))
  joined$stage_beta <- stage[, 1]
  joined$stage_p <- stage[, 2]
  joined$grade_beta <- grade[, 1]
  joined$grade_p <- grade[, 2]
  joined <- joined |>
    mutate(
      strict_candidate = main_fdr < THRESHOLDS$strict_survival_fdr &
        abs(main_beta) >= THRESHOLDS$min_abs_log_hr &
        abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc &
        abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc &
        sign(stage_beta) == sign(main_beta) & sign(grade_beta) == sign(main_beta) &
        stage_p < 0.05 & grade_p < 0.05,
      high_confidence = strict_candidate &
        main_fdr < THRESHOLDS$high_confidence_survival_fdr &
        abs(main_beta) >= THRESHOLDS$high_confidence_abs_log_hr
    )
  genes <- joined |> filter(high_confidence) |> pull(symbol)
  write_csv_atomic(joined |> filter(high_confidence) |> select(symbol, gene_id, main_beta, main_fdr),
                   file.path(out_dir, paste0(rule, "_high_confidence.csv")))
  tibble(rule = rule, n_patients = n_distinct(tumors_only$patient_barcode),
         n_reproducible = nrow(joined), n_strict = sum(joined$strict_candidate, na.rm = TRUE),
         n_high_confidence = length(genes),
         retained = paste(sort(intersect(genes, primary)), collapse = ";"),
         gained = paste(sort(setdiff(genes, primary)), collapse = ";"),
         lost = paste(sort(setdiff(primary, genes)), collapse = ";"))
}

summary_rows <- bind_rows(lapply(rules, fit_rule))
write_csv_atomic(summary_rows, file.path(DIRS$tables, "aliquot_sensitivity_summary.csv"))
message("Aliquot sensitivity complete.")
