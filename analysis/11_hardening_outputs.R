source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
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

writeLines(dossier_lines, file.path(DIRS$tables, "high_confidence_gene_dossiers.md"))

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
