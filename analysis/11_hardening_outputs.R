source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(survival)
  library(broom)
})

candidates <- read_csv(file.path(DIRS$tables, "candidate_gene_evidence_table.csv"), show_col_types = FALSE)
summary_counts <- read_csv(file.path(DIRS$tables, "candidate_summary.csv"), show_col_types = FALSE)

safe_neglog10 <- function(x) {
  out <- -log10(pmax(x, .Machine$double.xmin, na.rm = TRUE))
  out[!is.finite(out)] <- NA_real_
  out
}

rescale01 <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(1, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

gene_descriptions <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(candidates$symbol),
  keytype = "SYMBOL",
  columns = c("GENENAME", "ENTREZID")
) |>
  filter(!is.na(SYMBOL)) |>
  distinct(SYMBOL, .keep_all = TRUE) |>
  transmute(symbol = SYMBOL, gene_name = GENENAME, entrez_id = ENTREZID)

manual_context <- tribble(
  ~symbol, ~kidney_cell_type_sanity, ~literature_prior, ~manual_review_note,
  "CYFIP2", "Broad cytoskeletal/regulatory signal; not kidney-cell-type specific from this pipeline alone.", "ccRCC-specific role requires manual confirmation.", "Prioritize whether the survival signal reflects tumor biology or tissue composition.",
  "CRYL1", "Metabolic enzyme signal compatible with renal epithelial biology.", "Kidney cancer literature requires manual confirmation.", "Check whether downregulation is a renal identity/metabolic-loss signal.",
  "HHLA2", "Immune-checkpoint-like signal; may reflect tumor-immune interaction.", "Known cancer immunology gene; ccRCC-specific direction needs review.", "Interpret as immune/microenvironmental until cell source is resolved.",
  "DDC", "Metabolic/decarboxylase signal compatible with epithelial metabolic rewiring.", "Kidney cancer literature requires manual confirmation.", "Check whether expression is kidney-lineage, neuroendocrine-like, or stromal.",
  "KL", "Kidney-enriched aging/mineral-axis gene; plausible renal epithelial context.", "Strong kidney biology prior; ccRCC prognostic agreement needs review.", "A lower tumor expression signal may reflect loss of normal kidney program.",
  "PANK1", "Metabolism/coenzyme A pathway signal; plausible metabolic adaptation.", "Kidney cancer literature requires manual confirmation.", "Review whether signal is tumor-intrinsic metabolism or normal-tissue retention.",
  "IQGAP2", "Scaffold/cytoskeletal signaling; not cell-type-resolved here.", "Cancer literature exists; ccRCC-specific evidence needs review.", "Check direction against tumor suppressor versus composition explanations.",
  "ACADM", "Mitochondrial fatty-acid metabolism; plausible renal epithelial metabolism.", "Metabolic cancer relevance likely, ccRCC specificity needs review.", "Strong candidate for metabolic-adaptation interpretation.",
  "FUT6", "Glycosylation signal; cell type unresolved.", "Cancer glycosylation literature likely, ccRCC specificity needs review.", "Treat as pathway hypothesis until kidney/cancer literature is checked.",
  "CLCN5", "Renal proximal tubule/endocytic biology prior.", "Strong kidney biology prior; ccRCC prognostic agreement needs review.", "May mark retained normal tubule differentiation rather than aggressive tumor biology.",
  "C1QTNF6", "Secreted/metabolic-inflammatory signal; cell type unresolved.", "Cancer literature likely; ccRCC-specific role needs review.", "Check immune, adipokine, and tumor-cell sources.",
  "CADPS2", "Vesicle/secretion-related signal; not kidney-cell-type specific here.", "ccRCC literature requires manual confirmation.", "Lower-confidence biology unless supported externally.",
  "TEK", "Endothelial/angiogenesis-associated signal.", "Angiogenesis prior is relevant to ccRCC, but cell source is endothelial.", "Interpret as vascular composition or angiogenic state, not necessarily tumor-cell intrinsic.",
  "EMCN", "Endothelial marker-like signal.", "Kidney vascular literature likely; ccRCC survival interpretation needs review.", "Cell-type sanity flag: likely vascular composition.",
  "HIBCH", "Mitochondrial branched-chain amino-acid metabolism.", "Metabolic cancer relevance likely, ccRCC specificity needs review.", "Candidate metabolic-adaptation gene.",
  "TCIRG1", "Immune/lysosomal/proton-pump-associated signal; cell type unresolved.", "Immune/cancer literature likely, ccRCC specificity needs review.", "Check immune-cell composition before tumor-cell interpretation.",
  "PODXL", "Podocyte/endothelial/cell-adhesion signal.", "Cancer invasion literature likely; kidney compartment matters.", "Cell-type sanity flag: glomerular/endothelial biology may confound bulk signal.",
  "LRBA", "Immune-regulatory signal.", "Immune literature strong; ccRCC-specific direction needs review.", "Likely immune-composition sensitivity.",
  "GRAMD1A", "Cholesterol transport/contact-site biology; cell type unresolved.", "ccRCC literature requires manual confirmation.", "Review with lipid/metabolic adaptation lens.",
  "IFFO1", "Poorly resolved structural/nuclear signal.", "ccRCC literature likely sparse.", "Keep as computational association unless external evidence appears.",
  "ACAT1", "Mitochondrial acetyl-CoA/ketone and amino-acid metabolism.", "Metabolic cancer relevance likely, ccRCC specificity needs review.", "Strong candidate for metabolic-adaptation interpretation.",
  "DBT", "Branched-chain amino-acid metabolism.", "Metabolic cancer relevance likely, ccRCC specificity needs review.", "Strong candidate for metabolic-adaptation interpretation.",
  "TNFAIP2", "TNF/inflammatory-response signal.", "Inflammation/cancer literature likely; ccRCC specificity needs review.", "Interpret as immune/inflammatory until cell source is resolved.",
  "FHOD1", "Actin/cytoskeletal remodeling signal.", "Cancer migration literature likely; ccRCC specificity needs review.", "Check whether survival signal reflects invasion biology or stromal composition."
)

high_conf <- candidates |>
  filter(high_confidence_candidate) |>
  left_join(gene_descriptions, by = "symbol") |>
  left_join(manual_context, by = "symbol") |>
  mutate(
    reproducibility_component = rescale01(
      pmin(abs(tcga_log2fc), 5) / 5 +
        pmin(abs(gse40435_log2fc), 3) / 3 +
        pmin(abs(gse53757_log2fc), 3) / 3 +
        nominal_support_count / 2
    ),
    survival_strength_component = rescale01(safe_neglog10(main_fdr) + abs(main_log_hr)),
    ph_component = rescale01(pmin(main_ph_p_value, 1)),
    sensitivity_component = as.numeric(stage_sensitivity_same_direction & grade_sensitivity_same_direction) +
      as.numeric(stage_sensitivity_nominal & grade_sensitivity_nominal),
    pathway_component = if_else(pathway_class != "Unclassified", 1, 0),
    literature_component = case_when(
      grepl("Strong kidney biology prior|Angiogenesis prior|Immune literature strong", literature_prior) ~ 1,
      grepl("likely", literature_prior, ignore.case = TRUE) ~ 0.5,
      TRUE ~ 0
    ),
    final_rank_score =
      0.30 * reproducibility_component +
      0.30 * survival_strength_component +
      0.15 * ph_component +
      0.10 * sensitivity_component / 2 +
      0.10 * pathway_component +
      0.05 * literature_component,
    pubmed_ccrcc_query = paste0(
      "https://pubmed.ncbi.nlm.nih.gov/?term=",
      utils::URLencode(paste0(symbol, " (ccRCC OR clear cell renal cell carcinoma OR kidney cancer)"), reserved = TRUE)
    ),
    interpretation_status = "candidate prognostic association; requires manual literature and biology review"
  ) |>
  arrange(desc(final_rank_score), main_fdr)

ranked_shortlist <- high_conf |>
  transmute(
    rank = row_number(),
    symbol,
    gene_name,
    pathway_class,
    final_rank_score,
    evidence_score,
    tcga_log2fc,
    gse40435_log2fc,
    gse53757_log2fc,
    main_log_hr,
    main_hr,
    main_hr_ci_low,
    main_hr_ci_high,
    main_fdr,
    main_ph_p_value,
    stage_complete_log_hr,
    grade_complete_log_hr,
    kidney_cell_type_sanity,
    literature_prior,
    manual_review_note,
    pubmed_ccrcc_query,
    interpretation_status
  )

write_csv_atomic(ranked_shortlist, file.path(DIRS$tables, "high_confidence_ranked_shortlist.csv"))

manual_priority <- tribble(
  ~symbol, ~manual_tier, ~manuscript_role,
  "KL", "lead", "renal epithelial/metabolic retention",
  "ACADM", "lead", "fatty-acid oxidation and mitochondrial metabolism",
  "CRYL1", "lead", "renal metabolic differentiation",
  "ACAT1", "lead", "ketone, acetyl-CoA, and mitochondrial metabolism",
  "DDC", "lead", "amino-acid and biogenic-amine metabolism",
  "PANK1", "supporting", "CoA metabolism and mitochondrial function",
  "DBT", "supporting", "branched-chain amino-acid metabolism",
  "CLCN5", "supporting", "proximal-tubule endosomal/metabolic differentiation",
  "TCIRG1", "supporting risk", "lysosomal acidification, glycolytic/immune risk",
  "HHLA2", "interpret cautiously", "immune-checkpoint paradox",
  "GRAMD1A", "interpret cautiously", "cholesterol-contact-site biology",
  "C1QTNF6", "interpret cautiously", "secreted inflammatory/metabolic risk",
  "CYFIP2", "interpret cautiously", "cytoskeletal and broad tumor-suppressor-like signal",
  "IQGAP2", "interpret cautiously", "cytoskeletal scaffold signal",
  "TEK", "composition flag", "endothelial/vascular composition",
  "EMCN", "composition flag", "endothelial/vascular composition",
  "PODXL", "composition flag", "podocyte/endothelial/renal compartment signal",
  "FHOD1", "composition flag", "stromal/EMT and actin-remodeling signal",
  "IFFO1", "do not highlight", "weak ccRCC-specific biological support",
  "CADPS2", "do not highlight", "weak ccRCC-specific biological support",
  "LRBA", "do not highlight", "immune-composition signal",
  "FUT6", "do not highlight", "weak and directionally fragile glycosylation signal",
  "HIBCH", "do not highlight", "plausible but under-validated metabolic signal",
  "TNFAIP2", "do not highlight", "non-specific inflammatory signal"
)

manuscript_candidates <- ranked_shortlist |>
  left_join(manual_priority, by = "symbol") |>
  mutate(
    manual_tier = replace_na(manual_tier, "interpret cautiously"),
    manuscript_role = replace_na(manuscript_role, "candidate prognostic association requiring validation"),
    survival_direction = if_else(main_hr < 1, "higher expression associated with lower hazard", "higher expression associated with higher hazard"),
    tumor_direction = if_else(tcga_log2fc < 0, "lower in tumor", "higher in tumor")
  ) |>
  arrange(
    factor(manual_tier, levels = c("lead", "supporting", "supporting risk", "interpret cautiously", "composition flag", "do not highlight")),
    main_fdr
  )

write_csv_atomic(manuscript_candidates, file.path(DIRS$tables, "manuscript_candidate_prioritization.csv"))

fit_composition_sensitivity <- function(high_conf_symbols) {
  vst_mat <- read_required_rds(FILES$tcga_vst)
  coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
  clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
  repro <- read_csv(file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"), show_col_types = FALSE)

  ensembl <- sub("\\..*$", "", rownames(vst_mat))
  full_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ensembl),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  ) |>
    filter(!is.na(SYMBOL)) |>
    distinct(SYMBOL, .keep_all = TRUE) |>
    transmute(symbol = SYMBOL, gene_id = rownames(vst_mat)[match(ENSEMBL, ensembl)])

  clinical_surv <- clinical |>
    transmute(
      patient_barcode = submitter_id,
      os_time = os_time,
      os_event = os_event,
      age_years = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
      sex = gender,
      stage_clean = normalize_stage(ajcc_pathologic_stage),
      grade_clean = normalize_grade(tumor_grade)
    )

  tumor_samples <- coldata |>
    filter(sample_type == "Primary Tumor") |>
    mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
    inner_join(clinical_surv, by = "patient_barcode") |>
    filter(!is.na(os_time), os_time > 0, !is.na(os_event), sample_barcode %in% colnames(vst_mat))

  sample_barcodes <- tumor_samples$sample_barcode

  marker_sets <- list(
    proximal_tubule = c("AQP1", "LRP2", "CUBN", "SLC5A2", "ALDOB"),
    endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "TEK"),
    immune = c("PTPRC", "CD3D", "CD8A", "MS4A1", "CD68"),
    stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2")
  )

  score_df <- tibble(sample_barcode = sample_barcodes)
  marker_availability <- bind_rows(lapply(names(marker_sets), function(score_name) {
    present <- intersect(marker_sets[[score_name]], full_map$symbol)
    gene_ids <- full_map$gene_id[match(present, full_map$symbol)]
    gene_ids <- gene_ids[!is.na(gene_ids)]
    if (length(gene_ids) == 0) {
      score_df[[score_name]] <<- NA_real_
    } else {
      mat <- vst_mat[gene_ids, sample_barcodes, drop = FALSE]
      score_df[[score_name]] <<- as.numeric(scale(colMeans(t(scale(t(mat))), na.rm = TRUE)))
    }
    tibble(score = score_name, requested_markers = paste(marker_sets[[score_name]], collapse = ";"), used_markers = paste(present, collapse = ";"), n_used = length(gene_ids))
  }))

  base_dat <- tumor_samples |>
    transmute(
      sample_barcode,
      os_time,
      os_event,
      age = age_years,
      sex = factor(sex),
      stage = factor(stage_clean),
      grade = factor(case_when(
        grade_clean %in% c("G1", "G2") ~ "Low grade",
        grade_clean %in% c("G3", "G4") ~ "High grade",
        TRUE ~ NA_character_
      ))
    ) |>
    left_join(score_df, by = "sample_barcode")

  fit_one <- function(symbol) {
    gene_id <- repro$tcga_gene_id[match(symbol, repro$symbol)]
    if (is.na(gene_id) || !gene_id %in% rownames(vst_mat)) return(NULL)

    dat <- base_dat |>
      mutate(expr = as.numeric(scale(vst_mat[gene_id, sample_barcode]))) |>
      dplyr::select(os_time, os_event, expr, age, sex, stage, grade, proximal_tubule, endothelial, immune, stromal) |>
      filter(if_all(everything(), ~ !is.na(.x))) |>
      mutate(across(c(sex, stage, grade), droplevels))

    if (nrow(dat) < 100 || sum(dat$os_event == 1) < 25) return(NULL)
    if (any(vapply(dat[c("sex", "stage", "grade")], nlevels, integer(1)) < 2)) return(NULL)

    clinical_fit <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = dat)
    gene_fit <- coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade, data = dat)
    composition_fit <- coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade + proximal_tubule + endothelial + immune + stromal, data = dat)
    gene_term <- tidy(gene_fit, conf.int = TRUE) |> filter(term == "expr")
    comp_term <- tidy(composition_fit, conf.int = TRUE) |> filter(term == "expr")
    lrt_gene <- anova(clinical_fit, gene_fit, test = "LRT")

    tibble(
      symbol = symbol,
      n = gene_fit$n,
      events = gene_fit$nevent,
      clinical_concordance = unname(summary(clinical_fit)$concordance[1]),
      gene_model_concordance = unname(summary(gene_fit)$concordance[1]),
      composition_model_concordance = unname(summary(composition_fit)$concordance[1]),
      delta_concordance_gene_vs_clinical = gene_model_concordance - clinical_concordance,
      gene_lrt_p_vs_clinical = lrt_gene[2, "Pr(>|Chi|)"],
      base_log_hr = gene_term$estimate,
      base_hr = exp(gene_term$estimate),
      base_p_value = gene_term$p.value,
      composition_adjusted_log_hr = comp_term$estimate,
      composition_adjusted_hr = exp(comp_term$estimate),
      composition_adjusted_p_value = comp_term$p.value,
      same_direction_after_composition = sign(base_log_hr) == sign(composition_adjusted_log_hr)
    )
  }

  sensitivity <- bind_rows(lapply(high_conf_symbols, fit_one)) |>
    mutate(
      gene_lrt_fdr_vs_clinical = p.adjust(gene_lrt_p_vs_clinical, method = "BH"),
      composition_adjusted_fdr = p.adjust(composition_adjusted_p_value, method = "BH")
    ) |>
    arrange(gene_lrt_fdr_vs_clinical, composition_adjusted_fdr)

  write_csv_atomic(sensitivity, file.path(DIRS$tables, "candidate_clinical_composition_sensitivity.csv"))
  write_csv_atomic(marker_availability, file.path(DIRS$tables, "composition_marker_score_availability.csv"))
}

fit_composition_sensitivity(ranked_shortlist$symbol)

survival_report <- candidates |>
  filter(strict_candidate | high_confidence_candidate) |>
  transmute(
    symbol,
    tier = if_else(high_confidence_candidate, "high_confidence", "strict_only"),
    pathway_class,
    main_hr,
    main_hr_ci_low,
    main_hr_ci_high,
    main_log_hr,
    main_fdr,
    main_ph_p_value,
    main_n,
    main_events,
    ph_status = if_else(ph_pass, "pass", "fail"),
    main_warning_count,
    main_warning_text,
    stage_complete_log_hr,
    stage_complete_p_value,
    grade_complete_log_hr,
    grade_complete_p_value,
    sensitivity_pass
  ) |>
  arrange(desc(tier), main_fdr)

write_csv_atomic(survival_report, file.path(DIRS$tables, "candidate_survival_report.csv"))

threshold_grid <- tidyr::crossing(
  survival_fdr = c(0.01, 0.025, 0.05),
  min_abs_log_hr = c(log(1.25), log(1.5), log(2)),
  min_geo_abs_log2fc = c(0.25, 0.5)
)

threshold_sensitivity <- threshold_grid |>
  rowwise() |>
  mutate(
    n_candidates = sum(
      candidates$reproducible_deg &
        !is.na(candidates$main_fdr) &
        candidates$main_fdr < survival_fdr &
        candidates$ph_pass &
        abs(candidates$main_log_hr) >= min_abs_log_hr &
        abs(candidates$gse40435_log2fc) >= min_geo_abs_log2fc &
        abs(candidates$gse53757_log2fc) >= min_geo_abs_log2fc &
        candidates$sensitivity_pass,
      na.rm = TRUE
    )
  ) |>
  ungroup() |>
  arrange(survival_fdr, desc(min_abs_log_hr), desc(min_geo_abs_log2fc))

write_csv_atomic(threshold_sensitivity, file.path(DIRS$tables, "threshold_sensitivity.csv"))

set.seed(20260530)
eligible <- candidates |> filter(reproducible_deg, !is.na(main_fdr))
expression_pass <- eligible$geo_effect_support
survival_pass <- eligible$prognostic &
  eligible$ph_pass &
  eligible$meaningful_survival_effect &
  eligible$sensitivity_pass
observed_overlap <- sum(expression_pass & survival_pass, na.rm = TRUE)
n_perm <- 1000L
perm_overlap <- replicate(n_perm, {
  shuffled_survival_pass <- sample(survival_pass, length(survival_pass), replace = FALSE)
  sum(expression_pass & shuffled_survival_pass, na.rm = TRUE)
})

null_check <- tibble(
  check = "overlap_geo_effect_support_with_survival_filter_within_reproducible_deg",
  observed_overlap = observed_overlap,
  expected_overlap_mean = mean(perm_overlap),
  expected_overlap_sd = sd(perm_overlap),
  empirical_p_greater_equal = (sum(perm_overlap >= observed_overlap) + 1) / (n_perm + 1),
  n_permutations = n_perm,
  interpretation = "Permutation shuffles survival-filter pass labels within reproducible DEGs; it does not rerun Cox models."
)

write_csv_atomic(null_check, file.path(DIRS$tables, "null_overlap_check.csv"))

deg_rank <- candidates |>
  filter(reproducible_deg) |>
  arrange(desc(abs_tcga_log2fc)) |>
  slice_head(n = 25) |>
  transmute(symbol, list = "top_abs_tumor_normal_deg", rank = row_number(), abs_tcga_log2fc, abs_log_hr, main_fdr, strict_candidate, high_confidence_candidate)

prog_rank <- candidates |>
  filter(reproducible_deg, !is.na(abs_log_hr)) |>
  arrange(desc(abs_log_hr)) |>
  slice_head(n = 25) |>
  transmute(symbol, list = "top_abs_survival_effect", rank = row_number(), abs_tcga_log2fc, abs_log_hr, main_fdr, strict_candidate, high_confidence_candidate)

comparison <- bind_rows(deg_rank, prog_rank) |>
  group_by(symbol) |>
  mutate(appears_in_both_lists = n_distinct(list) > 1) |>
  ungroup() |>
  arrange(list, rank)

write_csv_atomic(comparison, file.path(DIRS$tables, "deg_vs_prognostic_comparison.csv"))

cell_type_sanity <- high_conf |>
  transmute(
    symbol,
    pathway_class,
    kidney_cell_type_sanity,
    bulk_rna_seq_risk = case_when(
      grepl("Endothelial|vascular|immune|podocyte", kidney_cell_type_sanity, ignore.case = TRUE) ~ "high_cell_composition_risk",
      grepl("renal epithelial|proximal tubule|metabolic", kidney_cell_type_sanity, ignore.case = TRUE) ~ "moderate_contextual_risk",
      TRUE ~ "unresolved_cell_source"
    ),
    manual_review_note
  )

write_csv_atomic(cell_type_sanity, file.path(DIRS$tables, "cell_type_sanity_check.csv"))

literature_table <- high_conf |>
  transmute(
    symbol,
    gene_name,
    pathway_class,
    literature_prior,
    pubmed_ccrcc_query,
    manual_review_note,
    curation_status = "seeded_for_manual_review"
  )

write_csv_atomic(literature_table, file.path(DIRS$tables, "high_confidence_literature_table.csv"))

dossier_lines <- c(
  "# HYDRA-ccRCC High-Confidence Gene Dossiers",
  "",
  "These dossiers are generated from the reproducible pipeline. They are interpretation scaffolds, not final biological claims.",
  ""
)

for (i in seq_len(nrow(ranked_shortlist))) {
  row <- ranked_shortlist[i, ]
  dossier_lines <- c(
    dossier_lines,
    paste0("## ", row$rank, ". ", row$symbol),
    "",
    paste0("- Gene name: ", ifelse(is.na(row$gene_name), "not available", row$gene_name)),
    paste0("- Pathway class: ", row$pathway_class),
    paste0("- Final rank score: ", round(row$final_rank_score, 3)),
    paste0("- Tumor-normal signal: TCGA log2FC ", round(row$tcga_log2fc, 3),
           ", GSE40435 log2FC ", round(row$gse40435_log2fc, 3),
           ", GSE53757 log2FC ", round(row$gse53757_log2fc, 3)),
    paste0("- Survival signal: HR ", round(row$main_hr, 3),
           " (95% CI ", round(row$main_hr_ci_low, 3), "-",
           round(row$main_hr_ci_high, 3), "), FDR ",
           format(row$main_fdr, scientific = TRUE, digits = 3),
           ", PH p ", round(row$main_ph_p_value, 3)),
    paste0("- Cell-type sanity: ", row$kidney_cell_type_sanity),
    paste0("- Literature prior: ", row$literature_prior),
    paste0("- Manual review note: ", row$manual_review_note),
    paste0("- PubMed query: ", row$pubmed_ccrcc_query),
    ""
  )
}

writeLines(dossier_lines, file.path("results", "high_confidence_gene_dossiers.md"))

funnel <- tibble(
  step = factor(
    c(
      "TCGA significant after symbol mapping",
      "Reproducible DEG",
      "Main Cox FDR < 0.05",
      "PH pass",
      "Sensitivity pass",
      "Strict candidate",
      "Pathway-classified high confidence",
      "High-confidence candidate"
    ),
    levels = c(
      "TCGA significant after symbol mapping",
      "Reproducible DEG",
      "Main Cox FDR < 0.05",
      "PH pass",
      "Sensitivity pass",
      "Strict candidate",
      "Pathway-classified high confidence",
      "High-confidence candidate"
    )
  ),
  count = c(
    sum(candidates$tcga_significant, na.rm = TRUE),
    summary_counts$value[match("reproducible_deg", summary_counts$metric)],
    summary_counts$value[match("main_stage_grade_complete_prognostic", summary_counts$metric)],
    summary_counts$value[match("ph_pass", summary_counts$metric)],
    summary_counts$value[match("sensitivity_pass", summary_counts$metric)],
    summary_counts$value[match("strict_candidate", summary_counts$metric)],
    sum(candidates$high_confidence_candidate & candidates$pathway_class != "Unclassified", na.rm = TRUE),
    summary_counts$value[match("high_confidence_candidate", summary_counts$metric)]
  )
)

write_csv_atomic(funnel, file.path(DIRS$tables, "evidence_funnel.csv"))

p_funnel <- ggplot(funnel, aes(step, count)) +
  geom_col(fill = "#285C7A", width = 0.68) +
  geom_text(aes(label = count), hjust = -0.12, size = 3.2) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "HYDRA-ccRCC Evidence Funnel",
    x = NULL,
    y = "Genes"
  ) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "evidence_funnel.png"), p_funnel, width = 7, height = 4.8, dpi = 300)

message("Hardening outputs complete.")
