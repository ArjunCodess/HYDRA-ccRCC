# Methods lock

Definitions below are the ones implemented in the analysis scripts. Sensitivities do not replace the primary rules.

## Cohort arithmetic

541 tumor aliquots and 533 patients differ by 8 aliquots. Those 8 are the extra samples from four patients with three tumor aliquots each. Normal aliquots are already one per patient: 72 normals, 72 patients, no patient with more than one normal. The 72 paired DESeq2 patients are those with both a selected tumor and a selected normal. The paired model is `~ patient + condition`.

## Sample choice

Within each patient and sample type the primary rule keeps the aliquot with the largest raw-count library depth. Ties break by barcode. That tie break is a deterministic reproducibility rule, not a biological claim. Depth was chosen so the retained file is the largest count library. `analysis/36_aliquot_sensitivity.R` refits the high-confidence gate after dropping the four multi-aliquot patients, after keeping the first barcode, and after summing a patient's tumor counts.

## Differential expression

A TCGA gene passes when FDR < 0.05 and the absolute MLE log2 fold change is at least 1. Nominal GEO support is a two-sided unadjusted p < 0.05 with the same sign as TCGA. A reproducible DEG is TCGA-significant, the same sign in both GEO cohorts, and nominal support in at least one GEO cohort.

GSE53757 has no sample-metadata field that gives the same patient id to every tumor and its normal. Adjacent rows are opposite tissues and the same stage in all 72 pairs. Four pairs also share a title token. The primary limma model remains `~ patient + condition`, with patient ids taken from that adjacent order. The unpaired model `~ condition` is a sensitivity. Its SVA null model is an intercept.

## Survival gates

The primary Cox model is overall survival on continuous expression standardized in the analyzed cohort, plus age, sex, AJCC stage, and collapsed grade. Stage is a categorical factor. With the default level order the reference is Stage I. Grade is G1+G2 versus G3+G4. In the TCGA clinical table, G1 contains 14 patients, G2 230, G3 207, and G4 78, so G1 is too small to keep as its own level. Age is linear. Patients missing age, sex, stage, or grade are excluded from the complete-case fits.

The stage sensitivity is `expr + age + sex + stage`. The grade sensitivity is `expr + age + sex + grade`. Both still adjust for age and sex.

Strict candidates also require main-model FDR < 0.05, absolute log hazard ratio at least log(1.25), absolute GEO log2 fold change at least 0.25 in each cohort, and nominal same-direction support in both sensitivity fits. High-confidence candidates also require FDR < 0.01 and absolute log hazard ratio at least log(1.5). These cutoffs are prespecified project thresholds in `analysis/00_config.R`. They are not a cited minimal clinically important difference.

The ridge arm is a training-fold glmnet Cox model. Clinical covariates are unpenalized. Every training-fold reproducible gene with a nonzero training standard deviation gets an L2 penalty (`alpha = 0`), so the fit does not select a sparse subset. `lambda.min` is chosen by partial-likelihood deviance on five event-stratified inner folds of the training patients. glmnet standardizes columns using the training matrix and applies those moments to the test rows. Expression is log2(count / size factor + 1). Training size factors come from the training DESeq2 fit. Test size factors use the training geometric-mean reference. The inner-fold seed is `RESAMPLING$seed + 31000 + repeat_id * 10 + fold`. Test outcomes are not an input. This arm is not the acceptance test.

The top gene is the highest evidence score:

`-log10(FDR) + |log HR| + min(|TCGA log2FC|, 5) / 5 + min(|GSE40435 log2FC|, 3) / 3 + min(|GSE53757 log2FC|, 3) / 3`

Ties break by FDR, then gene id. The prespecified prediction test adds that one gene. The ridge model on the training-fold reproducible set is a separate multigene description and is not the acceptance test.

## Nested validation

Ten repeats of five folds. Fold labels are redrawn at the start of each repeat with `set.seed(RESAMPLING$seed + 22 + repeat)` in both `analysis/22_nested_cv.R` and `analysis/31_nested_selection_benchmark.R`. Event-stratified means the labels are shuffled separately among deaths and among censored patients. Every outer test fold in the saved predictions contains 34 deaths.

Inside each training fold the procedure refits TCGA differential expression, applies the fixed GEO tables, refits the Cox gates, and chooses one gene. Library-size normalization and scaling for that gene are training-only. The GEO cohorts are not re-split. The 0.01 concordance rule is not re-chosen inside folds.

The patient-bootstrap interval resamples saved out-of-fold predictions. It does not refit selection or redraw folds. Call it a patient-bootstrap 95% interval conditional on the saved out-of-fold predictions.

Success requires a mean concordance increment of at least 0.01, a conditional interval above zero, and a mean three-year IPCW Brier difference less than or equal to zero. Exact zero on the Brier difference passes. The 0.01 increment is the prespecified size requirement for this project, not a universal clinical standard.

## Nulls

The exponential null draws event times from a clinical Cox cumulative hazard and reruns selection. It is a parametric clinical-risk null, not a p-value for the concordance increment. The permutation null shuffles training `(time, event)` pairs together, keeps expression fixed, reuses the five repeat-1 differential-expression fits, and cycles across those folds for 200 simulations. Seed `RESAMPLING$seed + 33`. The two frequencies are reported separately and are not pooled.

## Testing families

False-discovery correction is within one family at a time: TCGA differential-expression genes; each GEO cohort; main Cox models among reproducible DEGs; and, separately, the 23 high-confidence tests for shape, for composition, and for CheckMate interactions within each endpoint and model.

## Language

Prognostic means an association with outcome. Predictive is reserved for a treatment interaction. Tumor-normal differential expression is a disease-versus-kidney filter. It is not itself a prognostic discovery step. Their intersection can drop a prognostic gene that is not strongly differentially expressed. The all-gene Cox sensitivity, 12,630 associations at FDR < 0.05 among 32,192 genes, is that contrast.

## Marker scores

Each score is the mean of sample-standardized marker expression, then standardized across samples. Proximal-tubule markers are AQP1, LRP2, CUBN, SLC5A2, and ALDOB. Endothelial markers are PECAM1, VWF, KDR, EMCN, and TEK. Immune markers are PTPRC, CD3D, CD8A, MS4A1, and CD68. Stromal markers are COL1A1, COL1A2, DCN, LUM, and ACTA2. A candidate that is itself a marker is omitted from that score. These are bulk marker scores, not cell proportions. Consensus purity is an inferred bulk estimate. The six candidates that lose FDR support after marker adjustment are CLCN5, DDC, GJB1, HIBCH, KL, and PODXL. Purity adjustment does not remove those six.

## Gene mapping

| Dataset | Rule |
| --- | --- |
| TCGA-KIRC | Strip the Ensembl version, map with org.Hs.eg.db, keep the first symbol |
| GSE40435, GSE53757, GSE29609 | First token of the platform gene-symbol field, then the mean of probes that share it |
| E-MTAB-1980 | Strip the RefSeq version, map with org.Hs.eg.db, then the mean of rows that share a symbol |
| TRACERx Renal | Match the candidate symbol to expression rownames |
| CheckMate 025 | Mean of rows that share `gene_name` |

## TRACERx subset

The 39-patient subset is size-matched to GSE29609. It has 9 deaths. The reported percentages are medians across the 23 genes of each gene's direction-agreement rate. The interquartile ranges of those gene-level rates are part of the same summary. The histogram is the per-draw fraction of the 23 genes, which is a different summary of the same 1,000 draws.
