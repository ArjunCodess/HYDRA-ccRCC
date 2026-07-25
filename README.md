# HYDRA-ccRCC

**High-Discipline Reproducible Analysis of Clear Cell Renal Cell Carcinoma**

HYDRA-ccRCC is a reproducible public-data pipeline that tests whether genes consistently dysregulated between clear cell renal cell carcinoma (ccRCC) and normal kidney also carry stable, externally supported prognostic signal. The v2 pipeline uses TCGA-KIRC for discovery, GSE40435 and GSE53757 for tumor-normal replication, GSE29609 and E-MTAB-1980 for independent survival checks, repeated resampling for selection stability, repeated held-out testing for clinical increment, and Human Protein Atlas single-cell data for cell-source triangulation.

The completed run identified 3,304 reproducible differentially expressed genes, 519 strict TCGA-derived prognostic candidates, and a frozen 24-candidate high-confidence set. In E-MTAB-1980, 23 candidates mapped, 20 retained the TCGA hazard direction, and eight passed the prespecified external direction, FDR, adjusted-model, and proportional-hazards gates: `DDC`, `CRYL1`, `ACADM`, `TEK`, `EMCN`, `PODXL`, `FUT6`, and `HIBCH`.

This is an evidence-hardening study, not a validated biomarker panel. The strongest integrated metabolic leads are `ACADM`, `CRYL1`, and `DDC`; `TEK`, `EMCN`, and `PODXL` remain composition-sensitive vascular or renal-compartment signals; and `FUT6` and `HIBCH` need more biological validation despite external statistical support.

The development manuscript is in [`paper/main.pdf`](paper/main.pdf), with source in [`paper/main.tex`](paper/main.tex). It is an AI-assisted internal draft, not text for direct submission to IRIS or ISEF. The student must independently verify, rewrite, and cite the final submission.

## What v2 adds

- **Frozen independent survival validation:** E-MTAB-1980 tests all 24 candidates without feeding its results back into discovery or ranking.
- **Selection-stability analysis:** Twenty stratified 80% TCGA resamples rerun the complete 3,304-gene survival-selection screen.
- **Held-out clinical increment:** Twenty repeats of five-fold cross-validation compare clinical-only models with clinical-plus-one-gene models using out-of-fold concordance.
- **Cell-source triangulation:** Human Protein Atlas v25.1 normal-tissue single-cell data test whether candidate expression is vascular-, immune-, renal epithelial-, stromal-, or other-compartment dominant.
- **Reproducibility controls:** The pipeline records source URLs, access dates, file checksums, random seeds, package versions, and output-manifest checksums.
- **Automated validation:** Required columns, cohort sizes, funnel monotonicity, candidate coverage, and impossible values are checked before a run can finish.

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
| Proportional-hazards pass | 1,122 |
| Sensitivity-model pass | 1,063 |
| Strict candidates | 519 |
| Frozen high-confidence candidates | 24 |

### Independent survival checks

GSE29609 is a small microarray cohort with 39 samples and 17 deaths. Twenty-two of 24 candidates mapped, four retained the TCGA hazard direction, none had same-direction nominal support, and four had nominal opposite-direction associations (`DDC`, `ACAT1`, `CLCN5`, and `TCIRG1`). This cohort does not validate the shortlist, but its limited event count and platform coverage make it an exploratory contradiction check rather than a definitive refutation.

E-MTAB-1980 contains 101 patients and 23 deaths. Twenty-three candidates mapped, 20 retained the TCGA direction, 14 had same-direction unadjusted FDR support, 12 retained adjusted FDR support, and eight passed the complete external gate. No candidate had a nominally significant opposite-direction association.

The cohorts therefore give different answers. The larger E-MTAB-1980 result supports a subset of the frozen hypotheses, while the GSE29609 result prevents a claim of universal cross-platform replication.

### Stability and held-out performance

The stability analysis reran the full eligible-gene screen rather than refitting only the selected candidates. Six frozen candidates were selected in at least 80% of repeats: `PANK1` (100%), `ACADM` (95%), `DDC` (90%), `CLCN5` (85%), `CRYL1` (85%), and `CYFIP2` (85%). All 24 retained their full-fit hazard direction by median resampled effect, but the number selected per repeat ranged from 3 to 129, so the total threshold-defined set is unstable even when its strongest members are reproducible.

The repeated five-fold analysis produced a mean clinical-only concordance of 0.745. Adding each candidate individually gave a positive mean change for all 24; 23 had a positive empirical 2.5th percentile, while `FHOD1` crossed zero. The largest mean increments were for `CYFIP2`, `PANK1`, `CRYL1`, `ACADM`, `HHLA2`, and `CLCN5`. These estimates come from repeated splits of TCGA and do not establish clinical utility or validate a multigene model.

### Composition and interpretation

Bulk marker-score adjustment retained the TCGA hazard direction for all 24 candidates, but only 16 remained significant at FDR below 0.05. `KL` and `DDC` attenuated, which is consistent with part of their signal reflecting retained renal epithelial differentiation or tissue composition.

Human Protein Atlas normal-tissue single-cell data mapped all 24 candidates. `TEK` and `EMCN` were vascular-dominant; `CYFIP2`, `GRAMD1A`, and `TCIRG1` were immune-dominant under the prespecified mapping. These results reinforce composition cautions, but they do not replace tumor single-cell data or a direct tumor-purity estimate.

The combined v2 evidence changes the lead hierarchy:

- `ACADM`, `CRYL1`, and `DDC` have the strongest convergence of TCGA association, repeated selection stability, held-out increment, and E-MTAB-1980 support.
- `PANK1`, `CLCN5`, and `CYFIP2` are internally stable, but they do not pass the strict E-MTAB-1980 gate.
- `TEK`, `EMCN`, and `PODXL` pass the external gate but are best interpreted as vascular or renal-compartment signals.
- `FUT6` and `HIBCH` pass the external gate but remain biologically under-validated.
- `KL` and `ACAT1` remain biologically interesting TCGA-derived hypotheses, but neither receives strict E-MTAB-1980 support.

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

The latest complete run finished with:

```text
Pipeline complete.
```

`-SkipInstall` uses the existing local R library. `-ForceDownload` refreshes public inputs rather than using the cached processed files.

## Pipeline outputs

The pipeline generates:

- TCGA and GEO sample summaries and differential-expression tables
- cross-cohort reproducible DEG and evidence-funnel tables
- adjusted continuous-expression Cox models and sensitivity analyses
- high-confidence candidate tables, dossiers, and manual-review seeds
- GSE29609 and E-MTAB-1980 external survival results
- full-screen selection-stability results
- repeated held-out clinical-increment results
- bulk composition-sensitivity and HPA cell-source results
- GO enrichment, null checks, threshold sensitivity, and discordance analyses
- PCA, volcano, funnel, forest, and directional-discordance figures
- source provenance, package versions, and output checksums

Key v2 files:

- [`analysis/14_external_survival_emtab1980.R`](analysis/14_external_survival_emtab1980.R)
- [`analysis/15_selection_stability_tcga.R`](analysis/15_selection_stability_tcga.R)
- [`analysis/16_cv_clinical_increment.R`](analysis/16_cv_clinical_increment.R)
- [`analysis/17_hpa_cell_source.R`](analysis/17_hpa_cell_source.R)
- [`analysis/18_write_manifest.R`](analysis/18_write_manifest.R)
- [`analysis/12_validate_outputs.R`](analysis/12_validate_outputs.R)
- [`results/tables/external_survival_emtab1980.csv`](results/tables/external_survival_emtab1980.csv)
- [`results/tables/selection_stability_frozen_candidates.csv`](results/tables/selection_stability_frozen_candidates.csv)
- [`results/tables/candidate_cv_clinical_increment.csv`](results/tables/candidate_cv_clinical_increment.csv)
- [`results/tables/hpa_candidate_cell_source_summary.csv`](results/tables/hpa_candidate_cell_source_summary.csv)
- [`results/tables/source_provenance.csv`](results/tables/source_provenance.csv)
- [`results/tables/run_manifest.csv`](results/tables/run_manifest.csv)

## Dataset rules

All inputs are downloaded programmatically from public accession identifiers.

- TCGA-KIRC uses GDC STAR-count data through `TCGAbiolinks`.
- GSE40435 and GSE53757 provide independent tumor-normal expression replication.
- GSE29609 provides a small exploratory survival-direction check.
- E-MTAB-1980 provides the frozen primary external survival test.
- Human Protein Atlas v25.1 provides normal-tissue single-cell expression for cell-source triangulation.
- Raw FASTQ and BAM files are outside the project scope.
- GEO2R exports and manually edited metadata are not used.

## Scientific guardrails

- Continuous-expression Cox models are primary; median-split plots are secondary.
- The 24 candidates are frozen before E-MTAB-1980 testing.
- External null, missing, and contradictory results are retained.
- Repeated TCGA cross-validation estimates internal incremental discrimination, not clinical utility.
- Bulk marker scores are composition screens, not direct tumor-purity estimates.
- HPA data come from normal tissues, not ccRCC tumor single-cell samples.
- RNA associations do not establish protein effects, mechanisms, therapeutic targets, or causation.
- The strongest permitted label is **candidate prognostic association**.

## Completion state

The v2 computational pipeline and its prespecified analysis outputs are complete and pass automated validation. Direct tumor-purity estimation was not implemented, so composition conclusions remain limited to bulk marker scores and HPA triangulation.

Mentor review, student-authored IRIS materials, citation verification, source-data eligibility confirmation, and final submission are separate next-stage tasks. Internal mentor briefs and readiness checklists are intentionally kept out of version control.

AI/tool assistance is recorded in [`docs/AI_USE_LOG.md`](docs/AI_USE_LOG.md).
