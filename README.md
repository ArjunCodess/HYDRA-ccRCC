# HYDRA-ccRCC

HYDRA-ccRCC audits transcriptomic prognostic claims in clear cell renal cell carcinoma. The repository is an internal development study, not a validated clinical biomarker panel. [The manuscript](paper/main.pdf) is an AI-assisted draft that requires independent author verification before any submission.

The corrected primary analysis retains one TCGA primary tumor per patient: the aliquot with the highest raw-count library depth, with barcode as the tie-breaker. It includes 533 tumor patients and 72 selected normal patients; 72 patients have a selected tumor–normal pair. The previous 541-tumor-sample analysis and its 27-gene shortlist remain available in Git history. Regeneration produced 3,323 reproducible DEGs, 538 strict candidates, and 23 high-confidence candidates: [three genes entered, seven left, and 20 remained](results/tables/prior_candidate_delta.csv). This change matters because the old survival fits treated the three aliquots from each of four patients as independent observations.

The discovery cohorts are TCGA-KIRC, GSE40435, and GSE53757. All 23 candidates retained their tumor–normal direction in the 72-pair TCGA sensitivity, while 17 met its DE threshold. The already inspected GSE29609 and E-MTAB-1980 cohorts provide exploratory survival comparisons, not untouched external validation. E-MTAB-1980 mapped 22 of 23 candidates and gave strict same-direction support to 12; GSE29609 mapped 21 and agreed in direction for only five. DDC and TCIRG1 reversed direction with FDR below 0.05 in GSE29609, so a uniform cross-cohort replication claim is unsupported.

The matched-list funnel ablation compared the complete rule with survival-only, DE-only, both leave-one-GEO-out variants, and expression-matched controls at 22 genes per list. The complete rule did not outperform survival-only directional replication under the prespecified paired-interval criterion: the GSE29609 difference was 0.042 (95% bootstrap interval −0.194 to 0.270), and both rules achieved 100% direction agreement among their mapped E-MTAB-1980 genes. [The full gene results](results/tables/funnel_external_gene_results.csv) retain missing mappings, null results, and reversals.

<!-- nested-results:start -->
The selection-aware nested-CV result will be inserted from the completed result tables by `analysis/30_update_readme.R`.
<!-- nested-results:end -->

The corrected candidate set also has [conditional coefficient bootstraps](results/tables/candidate_cox_bootstrap_summary.csv), [conditional per-gene CV](results/tables/candidate_cv_clinical_increment.csv), matched-patient [purity](results/tables/candidate_direct_tumor_purity_sensitivity.csv) and [marker-score](results/tables/candidate_clinical_composition_sensitivity.csv) comparisons, [nonlinear and time-varying Cox checks](results/tables/candidate_survival_shape_sensitivity.csv), [TRACERx multiregion sensitivity](results/tables/tracerx_one_region_cox_summary.csv), and [CheckMate 025 treatment-interaction tests](results/tables/checkmate025_candidate_treatment_interactions.csv). The conditional bootstrap and per-gene CV reuse full-data selection, so they describe those 23 genes but cannot measure selection optimism. No CheckMate interaction survives multiplicity correction.

TRACERx contains 186 primary regions from 73 survival-linked patients, including 16 deaths. With a fixed event-stratified 39-patient subset, median candidate direction agreement across region draws was 98.5%; when patient membership also changed between repeats, it was 86.3%. The difference implicates patient-subset variation alongside region choice, and the nine events in each smaller subset limit precision. Human Protein Atlas data describe normal-tissue cell sources; they are not ccRCC tumor single-cell evidence.

## Reproduce the analysis

Run from the repository root on Windows with R 4.6.1 and MiKTeX:

```powershell
.\run_pipeline.ps1 -SkipInstall
.\build_paper.ps1
```

The local library's exact versions are recorded in [the package lock](environment/package_versions.csv). The cached public inputs have 638 paths, sizes, and hashes in [the input manifest](results/tables/input_manifest.csv); their original retrieval dates were not recorded. `analysis/00_verify_inputs.R` checks the cached hashes, and `analysis/12_validate_outputs.R` checks patient counts, pair coverage, candidate coverage, and FDR calculations. `-ForceDownload` refreshes public inputs and regenerates the manifest; it is not a reproduction from the committed cached inputs.

The full pipeline runs the 10-repeat, five-fold patient-level nested selection analysis and 200 clinical-only null simulations in [analysis/22_nested_cv.R](analysis/22_nested_cv.R). Each outer training fold repeats TCGA DE, GEO-supported filtering, adjusted survival selection, and training-only scaling. A fold without a selected gene falls back to the clinical-only prediction. Concordance, three-year IPCW Brier score, calibration slope, and a patient-level bootstrap interval are in the `nested_cv_*` tables. The [protocol](protocol.md) states the statistical rules and prespecified acceptance criteria.

The [run manifest](results/tables/run_manifest.csv) records generated-file checksums. All underlying public data stay under ignored `data/`; the PR commits code, derived tables and figures, source hashes, the protocol, and the development PDF. The strongest supported gene-level language is *candidate prognostic association*. New untouched external patients, tumor single-cell evidence, protein validation, and an independent randomized treatment cohort remain necessary for stronger claims.
