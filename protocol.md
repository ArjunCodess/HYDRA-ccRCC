# HYDRA-ccRCC Analysis Protocol

## Research question

Among genes reproducibly dysregulated across independent ccRCC expression cohorts, which have clinically adjusted survival associations that remain credible after sensitivity analysis, coefficient-uncertainty estimation, external outcome checks, and cell-source triangulation?

## Datasets

- TCGA-KIRC RNA-seq STAR counts and clinical data provide discovery expression and survival data.
- The primary TCGA table retains the highest raw-count library-depth tumor aliquot per patient, breaking ties by barcode; the same rule is applied to normal aliquots. The paired tumor--normal DE sensitivity uses patients with both selected sample types.
- GSE40435 and GSE53757 provide paired tumor-normal expression replication.
- GSE29609 provides a small external survival-direction stress test.
- E-MTAB-1980 provides the larger external survival evaluation.
- HPA v25.1 and Aran et al. consensus purity estimates support composition analyses.

The current candidate definition is a reviewer-driven reanalysis. External outcomes do not enter the selection rule, but the cohorts had been inspected in earlier project versions, so the revised set is not described as prospectively frozen or blindly validated.

## Differential expression and reproducibility

- TCGA uses DESeq2 on raw STAR unstranded counts. Significance requires FDR below 0.05 and absolute log2 fold change of at least 1.
- A separate sensitivity analysis estimates TCGA MAP log2 fold changes with `lfcShrink(type = "apeglm")`, applies no absolute fold-change inclusion threshold, and preserves the original MLE-based rule as the primary analysis.
- Each GEO cohort uses limma with patient-pair blocking. SVA protects the tumor-normal contrast using a full `patient + condition` model and a null `patient` model; estimated surrogate variables are added to the limma design.
- GSE53757 pair IDs are inferred from alternating row order. A condition-alternation assertion detects misordered rows, but it cannot independently establish patient identity without explicit subject identifiers.
- GEO tables report log2 fold-change confidence intervals and SVA design diagnostics. Zero estimated surrogate variables is retained as a valid result.
- A reproducible DEG must be TCGA-significant, have the same effect direction in both GEO cohorts, and have nominal p below 0.05 in at least one GEO cohort.

## Survival selection

- The outcome is overall survival, and expression is continuous and standardized.
- The main Cox model adjusts for age, sex, stage, and collapsed grade. Stage-complete and grade-complete models provide sensitivity checks.
- A strict candidate requires reproducible differential expression, main-model FDR below 0.05, absolute log hazard ratio of at least log(1.25), non-trivial GEO effects, and same-direction nominal support in both sensitivity models.
- A high-confidence candidate additionally requires main-model FDR below 0.01 and absolute log hazard ratio of at least log(1.5).
- `cox.zph` results are reported diagnostically. They do not exclude candidates, determine external support, or contribute to ranking; coefficients with diagnostic non-proportionality are interpreted as average hazard effects.
- The hard-threshold sensitivity fits the main age-, sex-, stage-, and grade-adjusted Cox model to every count-QC gene represented in the VST matrix. Benjamini--Hochberg correction is applied once across the complete successfully modeled universe, and the output is not used to redefine candidate membership.

## Coefficient uncertainty and prediction

- One thousand event-stratified patient bootstraps refit the age-, sex-, stage-, and grade-adjusted Cox coefficient for every high-confidence candidate.
- Bootstrap outputs include empirical standard error, percentile 95% interval, bias, direction agreement, and complete candidate-by-repeat results.
- Bootstrap p-values, per-repeat FDR thresholds, and selection frequencies are not used.
- Twenty repeats of five-fold event-stratified cross-validation compare clinical-only and clinical-plus-one-gene models on held-out patients. Cross-validation is used only to estimate prediction discrimination.
- The 20-by-five-fold per-gene analysis above is conditional on the full-data candidate set and is descriptive. A separate 10-by-five-fold patient-level outer CV reruns TCGA DE, both-GEO replication, and all candidate gates on the 517 patients with complete clinical covariates, with training-only normalization and expression scaling. The full-data stage-only and grade-only sensitivities can use additional patients with one missing covariate, so the nested fits use a narrower common analytic population. The nested fit chooses the highest evidence-score gene, uses clinical-only predictions if none qualifies, and reports the mean of per-repeat concordance, three-year IPCW Brier score, risk-score calibration slope, and a patient-resampled interval for the concordance increment. DESeq2 Cook outlier replacement follows the primary setting; any runtime failure triggers a logged no-replacement retry for that fold.
- Two hundred clinical-only null simulations repeat the fold-level selection on event times from an exponential baseline calibrated to the fitted clinical Cox cumulative hazard, with resampled censoring. One of the five fixed outer folds is assessed per simulation, cycling folds, so this null distribution estimates selection behavior under that specific clinical-only data-generating process rather than another 10-by-five-fold CV interval.
- The patient bootstrap resamples the saved out-of-fold predictions with patient identity linked across repeats. It does not refit selection or redraw folds, so its interval omits training-set and split uncertainty.

## External survival evaluation

- GSE29609 uses univariable continuous-expression Cox models because it contains 39 samples and 17 events.
- E-MTAB-1980 uses an unadjusted primary model and a limited secondary model adjusting for age, high T stage, and metastatic status.
- Revised strict external support requires TCGA-concordant direction and FDR below 0.05 in both E-MTAB-1980 models.
- External PH tests remain visible diagnostics but are not pass/fail gates.
- Changes in support status caused by removing a PH gate are identified explicitly and are not described as new data or stronger replication.
- Every revised candidate is reported, including missing, unsupported, and contradictory results.

## Multiregion transportability sensitivity

- Public TRACERx Renal primary-tumor TPM data are linked to patient, region, and overall-survival annotations from a pinned repository commit.
- For each high-confidence candidate, the high/low readout threshold is the median of patient-level regional medians. A multiregion tumor is discordant when at least one sampled region lies below the threshold and another lies at or above it.
- One thousand repeats select one primary-tumor region per patient and fit an unadjusted continuous-expression Cox model. A second scenario repeats the analysis in event-stratified 39-patient subsets to estimate the extra dispersion associated with a cohort the size of GSE29609.
- A third scenario fixes one event-stratified 39-patient subset before all 1,000 region draws, separating region-choice variation from changing which patients enter the model.
- This analysis measures regional readout and coefficient instability. It does not turn TRACERx into a blinded validation cohort or determine that spatial sampling caused the GSE29609 disagreement.

## Randomized treatment-interaction analysis

- Public supplementary data from Braun et al. link normalized pretreatment RNA expression to treatment and outcomes for the CheckMate 025 randomized comparison of nivolumab with everolimus.
- The analysis retains the 250 RNA-linked CheckMate 025 patients, standardizes each candidate across the combined trial subset, and fits a common Cox model containing treatment, continuous expression, and their interaction. Overall survival is the primary endpoint; progression-free survival is secondary.
- An unadjusted model preserves the randomized comparison. A limited sensitivity model adjusts for age, sex, and MSKCC risk group. Benjamini--Hochberg correction is applied across all high-confidence candidate interaction tests within each endpoint and model.
- Interaction estimates test differential treatment association, while treatment-arm-specific expression coefficients remain descriptive. This retrospective, post hoc analysis does not convert a prognostic candidate into a validated treatment-selection biomarker.

## Composition and interpretation

- Candidate models are tested with proximal-tubule, endothelial, immune, and stromal marker scores.
- Clinical, gene, marker-score, and direct-purity comparisons use each candidate's matched complete-case patients and the same expression scale. Nonlinear expression and expression-by-log-time alternatives are tested against the linear proportional-hazards model with false-discovery-rate control.
- Published consensus tumor purity is added directly to the clinical Cox model.
- HPA normal-tissue single-cell expression is mapped to renal epithelial, vascular, immune, stromal, and other compartments.
- Composition attenuation is evidence about interpretation, not a reason to discard a result.
- Integrated interpretation weighs external agreement, cohort contradiction, coefficient uncertainty, held-out concordance, and composition together rather than treating the external-support flag as a biological ranking.
- RNA-only associations are not described as mechanisms, therapeutic targets, or validated biomarkers.

## Reproducibility and acceptance

- Cached public inputs have committed checksums and package versions. Original retrieval dates were not recorded, so provenance does not substitute the pipeline run date for an access date.
- At matched list size, the full funnel is compared with survival-only, DE-only, each leave-one-GEO-out rule, and expression-matched control genes in the already inspected external cohorts. These are exploratory comparisons, not untouched validation.
- Engineering acceptance requires endpoint and patient-identity assertions, reproducible committed tables, and a built PDF. Scientific acceptance additionally requires a selection-aware CV concordance increment of at least 0.01 with a patient-level 95% interval above zero and no worse three-year Brier score; the external funnel comparison requires a paired 95% interval above zero against survival-only. Failed criteria are reported, not hidden.
- Automated validation checks cohort sizes, SVA diagnostics, candidate coverage, complete bootstrap grids, interval validity, output schemas, and funnel monotonicity.
- The strongest permitted gene-level label is **candidate prognostic association**.
