# HYDRA-ccRCC Analysis Protocol

## Research question

Among genes reproducibly dysregulated across independent ccRCC expression cohorts, which have clinically adjusted survival associations that remain credible after sensitivity analysis, coefficient-uncertainty estimation, external outcome checks, and cell-source triangulation?

## Datasets

- TCGA-KIRC RNA-seq STAR counts and clinical data provide discovery expression and survival data.
- GSE40435 and GSE53757 provide paired tumor-normal expression replication.
- GSE29609 provides a small external survival-direction stress test.
- E-MTAB-1980 provides the larger external survival evaluation.
- HPA v25.1 and Aran et al. consensus purity estimates support composition analyses.

The current candidate definition is a reviewer-driven reanalysis. External outcomes do not enter the selection rule, but the cohorts had been inspected in earlier project versions, so the revised set is not described as prospectively frozen or blindly validated.

## Differential expression and reproducibility

- TCGA uses DESeq2 on raw STAR unstranded counts. Significance requires FDR below 0.05 and absolute log2 fold change of at least 1.
- A separate sensitivity analysis estimates TCGA MAP log2 fold changes with `lfcShrink(type = "apeglm")`, applies no absolute fold-change inclusion threshold, and preserves the original MLE-based rule as the primary analysis.
- Each GEO cohort uses limma with patient-pair blocking. SVA protects the tumor-normal contrast using a full `patient + condition` model and a null `patient` model; estimated surrogate variables are added to the limma design.
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

## External survival evaluation

- GSE29609 uses univariable continuous-expression Cox models because it contains 39 samples and 17 events.
- E-MTAB-1980 uses an unadjusted primary model and a limited secondary model adjusting for age, high T stage, and metastatic status.
- Revised strict external support requires TCGA-concordant direction and FDR below 0.05 in both E-MTAB-1980 models.
- External PH tests remain visible diagnostics but are not pass/fail gates.
- Changes in support status caused by removing a PH gate are identified explicitly and are not described as new data or stronger replication.
- Every revised candidate is reported, including missing, unsupported, and contradictory results.

## Composition and interpretation

- Candidate models are tested with proximal-tubule, endothelial, immune, and stromal marker scores.
- Published consensus tumor purity is added directly to the clinical Cox model.
- HPA normal-tissue single-cell expression is mapped to renal epithelial, vascular, immune, stromal, and other compartments.
- Composition attenuation is evidence about interpretation, not a reason to discard a result.
- Integrated interpretation weighs external agreement, cohort contradiction, coefficient uncertainty, held-out concordance, and composition together rather than treating the external-support flag as a biological ranking.
- RNA-only associations are not described as mechanisms, therapeutic targets, or validated biomarkers.

## Reproducibility and acceptance

- Public sources, access dates, checksums, package versions, thresholds, and random seeds are recorded.
- Automated validation checks cohort sizes, SVA diagnostics, candidate coverage, complete bootstrap grids, interval validity, output schemas, and funnel monotonicity.
- The strongest permitted gene-level label is **candidate prognostic association**.
