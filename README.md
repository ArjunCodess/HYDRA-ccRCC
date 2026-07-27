# HYDRA-ccRCC

**High-Discipline Reproducible Analysis of Clear Cell Renal Cell Carcinoma**

HYDRA-ccRCC is a reproducible public-data pipeline that asks whether genes consistently dysregulated between clear cell renal cell carcinoma (ccRCC) and normal kidney also carry externally supported prognostic associations. It uses TCGA-KIRC for discovery, GSE40435 and GSE53757 for tumor-normal replication, GSE29609 and E-MTAB-1980 for survival checks, published consensus tumor-purity estimates, and Human Protein Atlas single-cell data for cell-source triangulation.

Following external statistical review, the analysis was revised to use surrogate-variable analysis (SVA) within the paired GEO models, report proportional-hazards tests as diagnostics rather than exclusion criteria, and use patient bootstrapping to estimate Cox coefficient uncertainty rather than repeated threshold-selection frequency. The resulting 27-candidate set is a reviewer-driven reanalysis. External outcomes were not used in its selection rule, but earlier versions of the project had already inspected those cohorts, so this repository does not describe the revised set as prospectively frozen or blindly validated.

The completed run identified 3,304 reproducible differentially expressed genes, 536 strict TCGA-derived prognostic candidates, and 27 high-confidence candidates. E-MTAB-1980 mapped 26 candidates and supported 13 under the revised same-direction, unadjusted-FDR, and adjusted-FDR rule: `DDC`, `CRYL1`, `ACADM`, `KL`, `ACAT1`, `CLCN5`, `TCIRG1`, `GJB1`, `TEK`, `EMCN`, `PODXL`, `FUT6`, and `HIBCH`.

This is an evidence-hardening study, not a validated biomarker panel. The development manuscript is in [`paper/main.pdf`](paper/main.pdf), with source in [`paper/main.tex`](paper/main.tex).

## What the pipeline includes

- **Paired, SVA-aware replication:** GSE40435 and GSE53757 use patient blocking, protect the tumor-normal contrast, estimate latent factors with SVA, and include any estimated surrogate variables in limma. The current run estimated zero additional surrogate variables in both cohorts.
- **Diagnostic PH testing:** `cox.zph` results remain in survival tables, but they do not exclude genes, determine external support, or contribute to candidate ranking.
- **Coefficient uncertainty:** One thousand event-stratified patient bootstraps refit each high-confidence adjusted Cox coefficient and report bootstrap standard errors, percentile 95% intervals, bias, and direction agreement.
- **Held-out prediction:** Twenty repeats of five-fold cross-validation compare clinical-only models with clinical-plus-one-gene models using out-of-fold concordance. Cross-validation is used only for prediction assessment.
- **External survival checks:** GSE29609 retains missing and contradictory results; E-MTAB-1980 reports unadjusted and limited-adjustment associations for every revised candidate.
- **Composition checks:** Published TCGA consensus purity, bulk marker scores, and HPA v25.1 normal-tissue single-cell expression expose tissue-composition explanations.
- **Reproducibility controls:** The pipeline records source URLs, access dates, checksums, random seeds, package versions, and output-manifest checksums, then validates required schemas and candidate coverage.

## Results

### Evidence funnel

| Metric | Value |
| --- | ---: |
| TCGA primary tumor / normal samples | 541 / 72 |
| GSE40435 tumor / normal samples | 101 / 101 |
| GSE53757 tumor / normal samples | 72 / 72 |
| TCGA-significant genes after symbol mapping | 8,852 |
| Reproducible DEGs | 3,304 |
| Main Cox FDR candidates | 1,176 |
| Sensitivity-model pass | 1,113 |
| Strict candidates | 536 |
| High-confidence candidates | 27 |

`RBM47`, `GJB1`, and `LTB4R` are the three high-confidence additions created by removing PH-based exclusion. Their PH diagnostic p-values are below 0.05, so their Cox coefficients are interpreted as average hazard effects and the non-proportionality signal remains visible in the result tables.

### External survival checks

GSE29609 contains 39 samples and 17 deaths. Twenty-five of 27 candidates mapped, six retained the TCGA direction, none had same-direction nominal support, and five had nominal opposite-direction associations: `DDC`, `ACAT1`, `CLCN5`, `TCIRG1`, and `RBM47`. Its small event count and different platform make it an exploratory contradiction check rather than definitive gene-level validation.

E-MTAB-1980 contains 101 patients and 23 deaths. Twenty-six candidates mapped, 23 retained the TCGA direction, 15 had same-direction unadjusted FDR support, and 13 retained adjusted FDR support and met the revised external rule. No candidate had a nominally significant opposite-direction association. PH diagnostics are reported for both external models but are not support gates.

### Bootstrap uncertainty and held-out prediction

All 27 candidates completed all 1,000 patient bootstraps. Every median bootstrap coefficient retained its full-fit direction, and every percentile 95% interval excluded zero; bootstrap standard errors ranged from 0.067 to 0.105. These estimates describe coefficient uncertainty and do not create another selection threshold.

The repeated five-fold analysis produced a mean clinical-only concordance of 0.745. Adding each candidate individually gave a positive mean change for all 27; 26 had a positive empirical 2.5th percentile, while `FHOD1` crossed zero. The largest mean increments were for `CYFIP2`, `PANK1`, `CRYL1`, `ACADM`, `HHLA2`, and `CLCN5`. These internal estimates do not establish clinical utility or validate a multigene model.

### Composition and interpretation

Bulk marker-score adjustment retained the TCGA hazard direction for all 27 candidates, while 18 remained significant at FDR below 0.05. Published consensus purity estimates matched 516 complete-case TCGA tumors with 170 deaths; all 27 candidates retained direction and FDR support after purity adjustment, with the largest relative attenuation, 8.9%, observed for `GRAMD1A`.

HPA normal-tissue single-cell data mapped all 27 candidates. `TEK` and `EMCN` were vascular-dominant, while `CYFIP2`, `GRAMD1A`, `TCIRG1`, and the newly retained `RBM47` were immune-dominant under the prespecified mapping. These observations reinforce composition cautions but do not replace tumor single-cell data.

The strongest integrated metabolic hypotheses remain `ACADM`, `CRYL1`, and `DDC`. Vascular, immune, inflammatory, and renal-compartment candidates require cell-source-aware interpretation even when their external statistics are supportive.

## Run the pipeline

From the repository root:

```powershell
.\run_pipeline.ps1
```

Useful options:

```powershell
.\run_pipeline.ps1 -SkipInstall
.\run_pipeline.ps1 -ForceDownload
```

`-SkipInstall` uses the existing local R library. `-ForceDownload` refreshes public inputs rather than using cached processed files.

## Key outputs

- [`analysis/05_deg_geo.R`](analysis/05_deg_geo.R) implements paired SVA-adjusted GEO differential expression.
- [`analysis/15_cox_bootstrap_uncertainty.R`](analysis/15_cox_bootstrap_uncertainty.R) estimates coefficient uncertainty without resample selection gates.
- [`analysis/16_cv_clinical_increment.R`](analysis/16_cv_clinical_increment.R) estimates held-out concordance.
- [`analysis/14_external_survival_emtab1980.R`](analysis/14_external_survival_emtab1980.R) performs revised external testing.
- [`results/tables/candidate_cox_bootstrap_summary.csv`](results/tables/candidate_cox_bootstrap_summary.csv) contains bootstrap uncertainty summaries.
- [`results/tables/candidate_cv_clinical_increment.csv`](results/tables/candidate_cv_clinical_increment.csv) contains cross-validated concordance changes.
- [`results/tables/external_survival_emtab1980.csv`](results/tables/external_survival_emtab1980.csv) retains external effect, FDR, and PH diagnostic results.
- [`results/tables/run_manifest.csv`](results/tables/run_manifest.csv) records generated-file checksums.

## Scientific guardrails

- Continuous-expression Cox models are primary; median splits are not used for selection.
- PH tests are diagnostics, and non-PH coefficients are interpreted as average hazard effects.
- External null, missing, and contradictory results are retained.
- Cross-validation estimates internal incremental discrimination, not clinical utility.
- Bootstrap intervals quantify coefficient uncertainty and do not define candidate membership.
- Published consensus purity is inferred rather than direct histopathologic measurement.
- HPA data come from normal tissues, not ccRCC tumor single-cell samples.
- RNA associations do not establish protein effects, mechanisms, therapeutic targets, or causation.
- The strongest permitted label is **candidate prognostic association**.

## Acknowledgment

We thank Levi Waldron for constructive statistical feedback on proportional-hazards screening, resampling, and batch adjustment. His review prompted the methodological reanalysis documented here; the authors remain solely responsible for the implementation, interpretation, and conclusions.

## Completion state

The revised computational pipeline completes from public accession identifiers and passes automated validation. Tumor single-cell, spatial, protein-level, prospective, and wet-lab validation remain outside the current pipeline.
