# HYDRA-ccRCC

**High-Discipline Reproducible Analysis of Clear Cell Renal Cell Carcinoma**

HYDRA-ccRCC tests whether genes dysregulated between clear cell renal cell carcinoma (ccRCC) and normal kidney also show prognostic associations. TCGA-KIRC supplies discovery data; GSE40435 and GSE53757 supply paired tumor-normal replication. Previously inspected GSE29609 and E-MTAB-1980 supply exploratory survival checks. Purity, marker scores, TRACERx Renal, and CheckMate 025 probe alternative explanations and transportability. This is an internal development study, not a validated biomarker panel. The [development manuscript](paper/main.pdf) is an AI-assisted draft that requires independent author verification before submission.

The corrected primary analysis selects one TCGA primary tumor per patient by highest raw-count library depth, breaking ties by barcode, and applies the same rule to normals. The earlier analysis counted three aliquots from each of four patients as independent tumors. Its 541 tumor samples and 27-gene set remain in Git history; the corrected analysis uses 533 tumor patients and yields 23 high-confidence candidates. [The candidate delta](results/tables/prior_candidate_delta.csv) records three additions, seven removals, and 20 retained genes. External outcomes do not enter the selection rule, but earlier versions had already inspected those cohorts, so they cannot provide untouched validation.

## What the pipeline includes

- **Patient-level discovery:** TCGA differential expression and survival analyses use one selected tumor per patient. A 72-pair TCGA tumor-normal analysis checks the discovery directions.
- **Paired GEO replication:** GSE40435 has 101 tumor-normal pairs and GSE53757 has 72. Limma models block on pairs and use surrogate variables if estimated; GSE53757 pair IDs are inferred from sample order rather than independently verified patient identifiers.
- **Survival sensitivity:** Proportional-hazards tests are diagnostics rather than selection gates. An all-gene `apeglm` analysis removes the hard fold-change gate and applies one global survival FDR correction.
- **Selection-aware prediction:** Ten repeats of five outer patient-level folds rerun candidate selection and training-only scaling. The comparator is a clinical-only model; a fold with no selected gene uses its clinical-only prediction.
- **Conditional uncertainty:** A separate 1,000-bootstrap coefficient analysis and per-gene held-out comparison characterize the full-data selected candidates. They do not estimate selection optimism.
- **Alternative explanations and transportability:** Matched-patient purity and marker-score fits, nonlinear and time-varying survival checks, TRACERx region draws, and CheckMate 025 treatment interactions test narrower hypotheses.
- **Reproducibility controls:** Cached input hashes, package versions, output checksums, patient joins, pair counts, survival-event coding, candidate coverage, and FDR calculations are checked. Original retrieval dates for cached inputs were not recorded.

## Results

### Evidence funnel

| Metric | Corrected value |
| --- | ---: |
| TCGA selected primary tumor / normal patients | 533 / 72 |
| TCGA selected tumor-normal pairs | 72 |
| GSE40435 / GSE53757 tumor-normal pairs | 101 / 72 |
| TCGA-significant genes after symbol mapping | 8,534 |
| Reproducible DEGs | 3,323 |
| Main Cox FDR candidates | 1,186 |
| Sensitivity-model pass | 1,117 |
| Strict candidates | 538 |
| High-confidence candidates | 23 |

The [funnel table](results/tables/candidate_summary.csv) reports each gate, and the [gene-level evidence table](results/tables/high_confidence_candidate_evidence.csv) lists the 23 genes alphabetically. All 23 retain their tumor-normal direction in the paired TCGA sensitivity, while 17 meet its differential-expression threshold. These counts describe selection within the discovery data; they are not external replication rates.

### Hard-threshold sensitivity

All 32,192 count-QC genes were fit in the age-, sex-, stage-, and grade-adjusted Cox sensitivity analysis. One global BH correction retained 12,630 at survival FDR below 0.05. This broad association burden makes the [all-gene analysis](results/tables/tcga_kirc_apeglm_all_gene_survival_summary.csv) a check on the fold-change gate, not a replacement shortlist or evidence that 12,630 genes are useful biomarkers.

### External survival checks

GSE29609 has 39 patients and 17 deaths. Of 23 candidates, 21 mapped, five kept the TCGA hazard direction, none had same-direction nominal support, and four had nominal opposite-direction associations. `DDC` and `TCIRG1` reversed with FDR below 0.05. The small event count limits individual estimates, but these [contradictory results](results/tables/external_survival_gse29609_summary.csv) rule out uniform external replication.

E-MTAB-1980 has 101 patients and 23 deaths. Of 23 candidates, 22 mapped, 21 kept the TCGA direction, 13 had same-direction FDR support, and 12 met the strict unadjusted and limited-adjustment rule. This [more favorable cohort](results/tables/external_survival_emtab1980_summary.csv) does not cancel the GSE29609 reversals. Both cohorts were inspected before the present rule was finalized and are exploratory.

At 22 genes per rule, the [matched-list funnel comparison](results/tables/funnel_external_gene_results.csv) tested the complete rule against survival-only, DE-only, both leave-one-GEO-out variants, and expression-matched controls. The complete-minus-survival-only directional-rate difference was 0.042 in GSE29609, with a paired gene-bootstrap 95% interval from -0.194 to 0.270, and zero in E-MTAB-1980. The prespecified external-funnel criterion failed. Different mapping denominators and correlation among genes further limit these descriptive intervals.

### Bootstrap uncertainty and held-out prediction

All 23 full-data selected candidates had conditional patient-bootstrap coefficient intervals excluding zero. The separate conditional per-gene cross-validation gave positive mean concordance increments for all 23. Both analyses reuse full-data selection, so neither shows what a new patient can expect from the entire selection procedure.

<!-- nested-results:start -->
The [selection-aware nested CV](results/tables/nested_cv_summary.csv) used 517 patients across ten repeats of five folds. It selected no gene in 0 folds; 2 folds required a DESeq2 no-replacement retry. The top-ranked gene varied across 10 genes, with the most frequent selected in 22 folds. The mean selected-gene minus clinical concordance was +0.0039 (patient-resampled 95% interval -0.0088 to +0.0171), and the three-year Brier-score difference was -0.0009. The prespecified CV criterion **failed**. Mean risk-score calibration slopes were 0.86 for clinical-only and 0.76 for the selected-gene models, so discrimination and calibration did not move together. The [200 clinical-only null simulations](results/tables/nested_cv_clinical_null_summary.csv) selected a gene in 3 single-fold runs under their simulated clinical-risk model. The patient bootstrap keeps the fitted fold models fixed, so its interval omits training-set and split uncertainty.
<!-- nested-results:end -->

The criterion required a concordance increment of at least 0.01, a patient-level 95% interval above zero, and no worse three-year Brier score. The observed increment misses the size threshold and its interval crosses zero. The null simulations characterize gene selection under their specified clinical-risk model; they are not a p-value for the repeated-CV increment.

### Composition and interpretation

In 516 matched complete-case TCGA patients with 170 deaths, direct consensus-purity adjustment retained all 23 candidate hazard directions and FDR associations. Marker-score adjustment on matched patients removed FDR support for six candidates. Neither analysis resolves whether the bulk signal is tumor-intrinsic, because purity estimates and marker scores are proxies for composition. The [interpretation table](results/tables/candidate_interpretation_context.csv) records gene-specific caveats without ranking “lead” biomarkers.

No candidate had an FDR-significant nonlinear or time-varying survival diagnostic in the [shape sensitivity](results/tables/candidate_survival_shape_sensitivity.csv). That is limited evidence against these particular departures, not proof of proportional hazards. Human Protein Atlas cell-source information comes from normal tissue and cannot substitute for ccRCC tumor single-cell evidence. The [earlier 24-gene biological review](results/archive/legacy_biological_review_24_candidates.md) is retained as a superseded historical record.

### Multiregion transportability

TRACERx linked 186 primary-tumor regions to 73 survival-linked patients, including 16 deaths; 58 patients had multiple regions, and all 23 candidates mapped. Across 1,000 draws from a fixed event-stratified 39-patient subset with nine deaths, median candidate direction agreement was 98.5% when region choice varied. It was 86.3% when patient membership varied as well. The [one-region analysis](results/tables/tracerx_one_region_cox_summary.csv) therefore implicates cohort composition alongside region choice; it does not establish that regional sampling caused the GSE29609 reversals.

### Randomized treatment interaction

The CheckMate 025 RNA-linked subset contains 250 patients, 120 assigned nivolumab and 130 everolimus, with 191 overall-survival and 222 progression-free-survival events. All 23 candidates mapped. Three overall-survival interactions were nominal in the unadjusted and adjusted models, but none survived within-endpoint, within-model FDR correction for either endpoint. The [interaction results](results/tables/checkmate025_candidate_treatment_interactions.csv) do not support a treatment-selection claim; this previously inspected subset is not an independent randomized confirmation cohort.

## Run the pipeline

From the repository root on Windows with R 4.6.1 and MiKTeX:

```powershell
.\run_pipeline.ps1 -SkipInstall
.\build_paper.ps1
```

`-SkipInstall` uses the existing local R library. The exact recorded versions are in the [package lock](environment/package_versions.csv). Allow several gigabytes of free disk space because regeneration writes a temporary DESeq2 RDS before replacing its cache. `-ForceDownload` refreshes public inputs and the input manifest; it is not a reproduction from the cached inputs. The [protocol](protocol.md) specifies model and acceptance rules.

## Key outputs

- [Input manifest](results/tables/input_manifest.csv): paths, sizes, and hashes for 638 cached public inputs. [Input verification](analysis/00_verify_inputs.R) checks them before analysis.
- [Funnel and candidate evidence](results/tables/candidate_summary.csv): corrected selection counts, with gene-level associations in the [candidate table](results/tables/high_confidence_candidate_evidence.csv).
- [Nested selection code](analysis/22_nested_cv.R) and [prediction results](results/tables/nested_cv_summary.csv): patient-level outer CV, clinical comparison, Brier scores, calibration, and uncertainty.
- [Funnel ablation results](results/tables/funnel_external_gene_results.csv): every rule, external mapping, null result, and reversal.
- [Conditional bootstrap](results/tables/candidate_cox_bootstrap_summary.csv) and [per-gene CV](results/tables/candidate_cv_clinical_increment.csv): descriptive analyses that retain the full-data selection condition.
- [Purity](results/tables/candidate_direct_tumor_purity_sensitivity.csv), [marker scores](results/tables/candidate_clinical_composition_sensitivity.csv), [TRACERx](results/tables/tracerx_candidate_multiregion_summary.csv), and [CheckMate 025](results/tables/checkmate025_candidate_treatment_interactions.csv): sensitivity and exploratory comparisons.
- [Acceptance criteria](results/tables/acceptance_criteria.csv) and [run manifest](results/tables/run_manifest.csv): scientific pass/fail status and generated-file checksums. [Output validation](analysis/12_validate_outputs.R) checks counts, joins, candidate coverage, and FDR calculations.
- [Manuscript source](paper/main.tex) and [development PDF](paper/main.pdf): interpretations tied to the regenerated tables.

## Scientific guardrails

- Continuous-expression Cox models are primary; median splits do not select candidates.
- PH diagnostics do not gate candidates, and a nonsignificant diagnostic does not prove proportional hazards.
- Missing external mappings, null results, and opposite-direction associations remain visible.
- Nested CV measures the complete training-set selection procedure; the conditional bootstrap and per-gene CV do not.
- The all-gene `apeglm` analysis is a sensitivity, not a retrospectively expanded candidate set.
- The funnel and nested-CV scientific acceptance criteria both failed; neither is described as confirmed superiority.
- Consensus purity is an estimate, marker scores are proxies, and HPA normal-tissue cells are not tumor single-cell validation.
- RNA associations establish neither protein effects nor mechanisms, therapeutic targets, causation, or clinical utility.
- The strongest supported gene-level label is **candidate prognostic association**.

## Acknowledgment

The earlier development draft acknowledged Levi Waldron, Michael Love, Philip Saylor, Samra Turajlić, and David McDermott for feedback. The author must verify those attributions and the manuscript independently before submission; responsibility for the implementation and interpretation remains with the author.

## Completion state

The corrected cached-input pipeline, output checks, and development PDF build complete with the recorded R 4.6.1 environment. The input manifest pins the cached files, but their original retrieval dates are unavailable. The prespecified prediction and external-funnel criteria failed. Untouched external patients, ccRCC tumor single-cell and protein evidence, and an independent randomized treatment cohort remain unavailable for stronger validation claims.
