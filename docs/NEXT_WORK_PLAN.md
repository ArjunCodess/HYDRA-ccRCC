# What to do next

The long pre-submission review mixes three kinds of advice. Most of it should not be the next week of work.

Already done, so do not rebuild it: patient-level correction (533 tumors, 23 high-confidence genes), the 10×5 nested CV, the six-arm nested benchmark (clinical, HYDRA gene, survival-only, DE-only, ridge, matched control), the exponential clinical-only null, TRACERx region-versus-patient resampling, CheckMate 025 interactions, purity versus marker-score adjustment, and the matched-list external funnel comparison. The manuscript already reports the ridge concordance increment of about 0.026 with a Brier interval that still crosses zero. That was the review’s “single biggest experiment.”

Do not spend the week rewriting `paper/main.tex`. IRIS/ISEF requires the student to write the synopsis, paper, abstract, poster, and citations. The development draft stays an internal numbers source. Title, author line, “development manuscript” banners, and talking points for judges are a submission step after the results below are frozen.

Skip these. They do not change what the project can claim:

- Another penalized model (lasso or elastic net). Ridge already showed that the reproducible set carries held-out ranking information and that the one-gene rule does not.
- An oracle “best gene on the full dataset” run. The conditional per-gene CV already shows optimistic increments for genes chosen on all of TCGA. The nested increment of +0.0039 is the honest contrast.
- `renv.lock` and a new `make all`. `environment/package_versions.csv`, the input manifest, and `run_pipeline.ps1` already exist.
- Prose substitutions (“failed” to “did not meet the criterion,” “support” glossaries, HYDRA definition, novelty paragraph). Put the locked definitions in the methods sheet below, then use them when you write the submission yourself.

Five pieces of work are worth doing. Do them in this order. The first two can change the 23-gene set. The next two harden the validation claim. The last one makes the results a reviewer can check without hunting.

Stop rule: these sensitivities do not silently replace the primary 23. If one of them moves the shortlist, write that down and decide before rerunning the nested CV. The nested CV and the benchmark are expensive and should stay tied to the current primary rule unless you explicitly change the protocol.

---

## 1. Find out whether GSE53757 pairing is real

This is the sharpest methodological hole. `analysis/05_deg_geo.R` assigns GSE53757 patient IDs by alternating row order (`rep(seq_len(n/2), each = 2)`). The only check is that adjacent rows have different conditions. The reproducible-DEG gate treats those pairs as patients. A reviewer can reject the 3,323-gene list on that sentence alone.

### What to do

1. Load the cached GSE53757 ExpressionSet and print every `pData` column for the first 12 samples: `title`, `geo_accession`, `source_name_ch1`, `characteristics_ch1`, and any column whose values repeat in pairs. GSE40435 already parses `title` with `patient\\s+([0-9]+)`. GSE53757 may have the same kind of token and the parser never looked.
2. If a real subject field exists, rebuild pair IDs from it, assert one tumor and one normal per subject, and diff the limma table against `results/tables/gse53757_limma_tumor_vs_normal.csv`. Record how many of the 23 candidates change GEO direction, nominal p, or the 0.25 absolute log2 fold-change gate.
3. If no subject field exists, do not invent one. Add an unpaired sensitivity in the same script, gated by a flag so the primary paired table stays the primary table:
   - Design `~ condition` instead of `~ patient + condition`.
   - SVA null model `~ 1` instead of `~ patient`.
   - Same gene-collapsing function.
4. Join the unpaired GSE53757 statistics to the existing GSE40435 table and the TCGA DEG table. Recompute only the reproducibility gate: TCGA FDR < 0.05, absolute log2 fold change ≥ 1, same direction in both GEO cohorts, nominal p < 0.05 in at least one. Write `results/tables/gse53757_unpaired_sensitivity.csv` with one row per current high-confidence gene: paired log2FC, unpaired log2FC, paired p, unpaired p, still reproducible yes/no.
5. Add one summary row: of 23, how many still pass. Also report the Jaccard overlap of the full reproducible-DEG symbol lists, not just the 23.

### How to read it

- A real subject field that matches the alternating order: say that in the methods and retire the “inferred from order” limitation.
- A real subject field that disagrees: the primary analysis has to move to the metadata pairs. That is a protocol change. Stop and rerun the funnel before any other new claim.
- No subject field, and all or nearly all of the 23 still pass unpaired: keep the paired analysis as primary and report the unpaired overlap as the robustness result.
- No subject field, and several of the 23 drop: the reproducible gate is sensitive to an unverified pairing. That becomes a main limitation, not a footnote. Do not rerun nested CV until you decide whether the primary gate stays paired.

### Files

- Edit `analysis/05_deg_geo.R` only for the pData audit and the unpaired branch.
- New table consumed later by `analysis/25_paper_numbers.R` if you add a macro. Do not add manuscript prose beyond one sentence and the macro.

---

## 2. Show that the aliquot rule is not the result

The arithmetic is already consistent and should be written once so nobody re-derives it wrong: 541 tumor aliquots, 533 patients, 8 extra aliquots, and those 8 are exactly the extra samples from four patients with three aliquots each (4 patients × 2 extras = 8). Selection is in `analysis/functions/patient_samples.R`: within patient and sample type, sort by descending raw-count library depth, then barcode, and keep row 1. The barcode sort is a reproducibility tie-break, not a biological claim. Say that. The open question is whether “deepest library” chooses a different shortlist than the obvious alternatives.

Only four patients have a real choice. That makes a full recount of all 32,192 Cox models unnecessary.

### What to do

1. From the TCGA colData and counts, write `results/tables/tcga_aliquot_selection_audit.csv`: one row per tumor aliquot that belongs to a patient with more than one tumor aliquot. Columns: patient, barcode, sample type, library depth, whether it was selected. The table should show four patients and twelve aliquots.
2. Three alternate matrices, each still one row per patient, rerun only through the high-confidence gate (DESeq2, the two GEO gates as they are, the main Cox model and the two sensitivity fits). Do not rerun nested CV.
   - Drop the four multi-aliquot patients entirely.
   - Keep the lexicographically first barcode instead of the deepest library.
   - Sum raw counts across a patient’s tumor aliquots, and keep the existing single normal.
3. Implementation: add a `rule` argument next to `select_tcga_patient_samples` rather than copying the DESeq2 script three times. Rules: `deepest` (current), `first_barcode`, `drop_multi`, `sum_counts`. `analysis/04_deg_tcga.R` and `analysis/07_survival_tcga.R` should be callable on an alternate coldata without overwriting `tcga_kirc_deseq2_tumor_vs_normal.csv` or `tcga_kirc_cox_models.csv`. Write alternate tables under `results/tables/aliquot_sensitivity/`.
4. Summary table `results/tables/aliquot_sensitivity_summary.csv`: rule, n patients, n reproducible DEGs, n strict, n high-confidence, symbols gained, symbols lost, symbols retained, versus the primary 23.

### How to read it

- Same 23, or a loss of one or two genes, under all three rules: the depth rule is not driving the shortlist. One paragraph plus the table closes the review point. State in the methods that depth was chosen so the selected library is the largest count file, that ties break by barcode for determinism, and that the three alternates are the sensitivity.
- A different 23 under `sum_counts` or `first_barcode`: report it. Keep deepest-library as the prespecified primary rule. Do not switch rules after seeing survival results.

Also record, in the same audit table, the normal side: how many normal aliquots existed before selection, how many patients had more than one normal, and that the 72 pairs are the patients who have both a selected tumor and a selected normal. The paired DESeq2 model is `~ patient + condition` in `analysis/04b_paired_deg_tcga.R`. Cite that formula. Do not make the reader infer it.

---

## 3. Put a published signature on the same held-out patients

The nested benchmark shows that the funnel’s one gene does not beat the clinical model, and that a ridge on the funnel’s own gene pool does. A reviewer will still ask how either compares with a signature the field already uses. E-MTAB-1980 already carries a `molecular_subtype` column restricted to ccA/ccB in `analysis/14_external_survival_emtab1980.R` and `analysis/24_funnel_ablations.R`. That label was assigned by the original study, not by HYDRA.

Do not retrain ccA/ccB or ClearCode34 on TCGA. That would be another selection procedure.

### What to do

1. E-MTAB-1980, using the existing processed expression and clinical table:
   - Clinical-only C-index, three-year Brier, calibration slope, with the covariates already used there (unadjusted primary, and the limited age / high-T / metastatic model).
   - The same metrics for ccB versus ccA as a two-level factor. Patients without a ccA/ccB label stay out of this comparison only, and the excluded count is reported.
   - The same metrics for the clinical model plus the single HYDRA gene that the full-data rule would have added, and for a fixed sum of the mapped high-confidence genes using the TCGA log-hazard signs. The signs are frozen from `high_confidence_candidate_evidence.csv` before looking at E-MTAB performance.
   - This cohort was previously inspected. Label the table exploratory. It answers whether a published subtype separates patients in this cohort more clearly than the funnel’s genes. It is not a new discovery set.
2. TCGA, on the existing 517 complete-case patients and the existing out-of-fold clinical predictions in `results/tables/nested_cv_predictions.csv`:
   - Implement ClearCode34 only if you transcribe the 34 gene symbols and the published good-prognosis versus poor-prognosis direction from Brooks et al., European Urology, 2014, without refitting those directions.
   - Score each patient as the mean of training-fold standardized expression with those published signs. Standardization mean and sd come from the training fold only, using the same fold IDs as `analysis/22_nested_cv.R` (`RESAMPLING$seed + repeat`).
   - Add that score, unpenalized, to the clinical Cox model inside each training fold. Score the held-out fold with the existing concordance, Brier, and calibration functions.
   - This arm does not select genes. It does not need DESeq2. It should be a short script, `analysis/32_published_signature_benchmark.R`, that reads the saved fold assignments and counts, not a fork of the full nested CV.
3. Write `results/tables/published_signature_summary.csv` with one row per arm and cohort: n, events, C-index, delta versus clinical, Brier delta, calibration slope. Add a patient-bootstrap interval with the same caveat already used everywhere: conditional on the saved predictions, not a refit of selection.

### How to read it

- If ClearCode34 or ccA/ccB beats the clinical model and the HYDRA gene does not, the paper’s claim gets sharper: the failure is specific to this funnel’s one-gene rule, in a setting where a published score can still move discrimination.
- If neither published score beats the clinical model by the prespecified 0.01 either, say that. Small external cohorts (23 deaths in E-MTAB-1980) can make every molecular score look similar. Do not declare the published signature a winner or a loser from a noisy C-index.
- If you cannot transcribe ClearCode34 from the paper without guessing a gene or a sign, ship the E-MTAB ccA/ccB comparison only. A wrong gene list is worse than a missing arm.

---

## 4. Add a permutation null that does not refit DESeq2

The current null (`NullSelected = 3` of 200 in `analysis/22_nested_cv.R`) draws exponential event times from the fitted clinical cumulative hazard and reruns selection on one fixed fold per simulation. Two limits are already stated and should stay stated: it is not a p-value for the +0.0039 increment, and the baseline hazard is parametric.

The second null should break the expression–outcome link and leave the expression selection pool alone. Differential expression does not use survival, so rerunning DESeq2 inside a survival null wastes the run and does not answer a new question.

### What to do

1. New script `analysis/33_survival_permutation_null.R`.
2. Reuse one repeat of the existing event-stratified fold IDs (repeat 1). For each of 200 permutations, cycle across the five folds the same way the current null does.
3. Inside the training fold, permute the pairs `(os_time, os_event)` across training patients. Keep each patient’s expression and covariates fixed. Do not permute inside a clinical stratum unless the unstratified permutation is somehow invalid; unstratified permutation is the sharper null because it also breaks the clinical–outcome link that the exponential null preserves. Document which one you ran. Prefer the unstratified permutation and keep the exponential null as the “clinical risk is real, expression is not” contrast. Together they answer two different questions. Do not replace the exponential null.
4. Rank and test genes with the same high-confidence rule used in `select_gene`: FDR within the training reproducible set, absolute log hazard ratio at least `log(1.5)`, and the stage-only and grade-only checks. The reproducible set for that fold can be recomputed from the saved training DE if it is cached; if it is not cached, refit DESeq2 once per fold (five fits total) and reuse those five objects across all 200 permutations.
5. Record, per simulation: fold id, whether any gene passed, the symbol that would have been selected, and its training FDR. Summary: selection frequency, and the distribution of the selected gene’s training absolute log hazard ratio when one was selected.
6. Seed: `RESAMPLING$seed + 33`. Write it into the script header and into `protocol.md` in the same edit that describes the null.

### How to read it

- Selection in a small fraction of permutations means the high-confidence gate is not trivially easy to pass when outcomes are scrambled. It still does not test the out-of-fold C-index.
- Selection in a large fraction means the gate is loose once the reproducible set is fixed. That would qualify the “3/200” sentence rather than reinforce it.
- Report both frequencies side by side. Do not pool them into one p-value.

---

## 5. Make the existing numbers checkable

These use tables that already exist. Do them even if you stop after item 2.

### Cross-cohort forest

`results/tables/external_survival_gse29609.csv` and `results/tables/external_survival_emtab1980.csv` already contain `external_log_hr`, confidence limits, and `main_log_hr` for each candidate. `analysis/09_figures_tcga.R` draws a TCGA-only forest (`candidate_forest_plot.png`). Add `analysis/34_external_forest.R` that plots, for each of the 23 symbols, the TCGA log hazard ratio and the two external log hazard ratios with intervals. Mark DDC and TCIRG1. Order genes by the TCGA estimate. One figure, `results/figures/external_log_hr_forest.png`. This is the picture that stops a reader from treating “21 of 22 same direction” as replication. Direction, a confidence interval, and an FDR are different claims; the figure should show the intervals, and the caption should say that GSE29609 has 17 deaths and E-MTAB-1980 has 23.

### Fold event counts

`nested_cv_folds.csv` has sizes but not events. `nested_cv_predictions.csv` has `os_event`. Summarize events per repeat and fold: min, median, max in the test fold and in the training fold. Write `results/tables/nested_cv_fold_events.csv`. The methods sentence is then numerical: event-stratified means the fold labels were assigned by shuffling within death and within censoring separately (`folds()` in `analysis/22_nested_cv.R`), and the test folds contained this range of deaths. The five folds are redrawn each repeat with `set.seed(RESAMPLING$seed + 22 + repeat)` in the benchmark and the corresponding seed in `analysis/22_nested_cv.R`. Confirm those two seeds match before writing the sentence. The review asked whether folds are fixed across repeats. They are not. They are regenerated each repeat. Say that, and say the null cycles one fold from the fixed repeat-1 partition.

### Methods lock

Add `docs/METHODS_LOCK.md`, copied from the code rather than from memory. The submission draft should follow this sheet. Include only definitions that are currently true:

| Choice | Locked definition |
| --- | --- |
| Nominal support | Two-sided Wald p < 0.05, unadjusted, and the same sign as TCGA |
| Reproducible DEG | TCGA FDR < 0.05 and \|log2FC\| ≥ 1, same sign in both GEO cohorts, nominal support in at least one GEO cohort |
| Grade | G1+G2 versus G3+G4, because the script collapses those pairs in `analysis/22_nested_cv.R` and the primary Cox script. State the raw G1–G4 counts in the sheet so the sparsity reason is a count, not a slogan. If G1 is not sparse, say the collapse was a prespecified grouping and keep it. Do not uncollapse after seeing coefficients. |
| Stage | Categorical AJCC stage after `normalize_stage`, not a linear score. Reference level is whatever `factor()` drops first; write that level down after printing `levels()`. |
| Standardization | Training patients only, inside each nested fold. Full-data models standardize on the analyzed complete-case cohort. |
| HR gates | \|log HR\| ≥ log(1.25) for strict and ≥ log(1.5) for high confidence. These are prespecified analysis thresholds in `analysis/00_config.R`. They were not taken from a cited clinical minimally important difference. Say that. Do not hunt for a citation that does not match. |
| GEO effect gate | Absolute log2 fold change ≥ 0.25 in each GEO cohort. Same status: prespecified project threshold. |
| Stage-only / grade-only | The gene plus stage, or the gene plus grade. Age and sex are not in those two fits. Confirm against `analysis/07_survival_tcga.R` and paste the formulas. |
| “No worse Brier” | Mean gene-minus-clinical three-year IPCW Brier ≤ 0. Exact zero passes. No extra tolerance. |
| ΔC ≥ 0.01 | Prespecified size requirement in `docs/ACCEPTANCE_CRITERIA.md` and `protocol.md`, not a universal clinical standard. The interval must also lie above zero. |
| What “selection-aware” includes | Inside each training fold: TCGA tumor-versus-normal DESeq2, GEO gate using the fixed full GEO tables, Cox gates, ranking, and the choice of one gene. Training-only library-size normalization and scaling. |
| What it does not include | The GEO cohorts are not re-split. The 0.01 rule was not re-chosen inside folds. The ridge arm was specified after the one-gene rule and is not the acceptance test. The patient-bootstrap interval does not refit selection or redraw folds. Call it a patient-bootstrap interval conditional on the saved out-of-fold predictions, in the abstract as well as the methods. |
| Testing families | One row per family: TCGA DE genes; GEO genes within each cohort; main Cox models among reproducible DEGs; the 23 high-confidence tests for shape, for composition, and for CheckMate interactions within endpoint and model. Do not imply one FDR across all of these. |
| Top gene | Highest evidence score already used by `select_gene`. Ties break by gene id inside that score. State the score formula from the function, not a paraphrase. |
| Why one gene | The prespecified acceptance test is one added gene, to avoid a second selection layer over the 23. The ridge arm is the multigene answer and it is already computed. Point to it instead of promising a new multigene model. |
| Predictive versus prognostic | Prognostic means an association with outcome. Predictive means a treatment interaction. Only the CheckMate section may use “predictive,” and only for the interaction test. The current manuscript is already careful. Keep it that way in anything you write. |
| Tumor-normal versus survival | Tumor-normal DE is a disease-versus-kidney filter. Survival is a separate filter. Their intersection can drop a prognostic gene that is not strongly differentially expressed. The all-gene result (32,192 tested, 12,630 at FDR < 0.05) is the illustration. The nested benchmark is the evidence about whether that intersection helped prediction. You do not need a new analysis to make this point. |

While you are in the docs, correct two stale files that still say 27 candidates: `docs/MENTOR_BRIEF.md` and `docs/ACCEPTANCE_CRITERIA.md`. Point them at 23 and at the nested result. A mentor reading the brief today will review the wrong project.

### Cohort table

One CSV, `results/tables/cohort_dictionary.csv`, one row per dataset already in the pipeline: TCGA-KIRC tumors, TCGA normals, TCGA pairs, GSE40435, GSE53757, GSE29609, E-MTAB-1980, TRACERx, CheckMate 025. Columns: accession, role in this project, patients, samples or regions, events, tissue, platform, previously inspected yes/no, script that builds it. Fill from the summary tables and `README.md`. This is the table the review asked for. It does not require new statistics.

---

## Order and time

1. Methods lock, cohort dictionary, fold-event table, external forest. About a day. No model fitting. Do this first so the later results have a place to land.
2. GSE53757 pData inspection. An hour if the metadata is already in the cached ExpressionSet. The unpaired refit is the long part.
3. Aliquot audit table, then the three alternate gates. The drop-four-patients run is the one to finish even if `sum_counts` is slow.
4. Permutation null. Five DESeq2 fits, then 200 Cox screens. Much smaller than the original nested CV.
5. Published-signature arm. E-MTAB ccA/ccB first. ClearCode34 only after the gene list is transcribed.

After each item, add the summary numbers through `analysis/25_paper_numbers.R` and one results sentence. Do not rewrite the introduction.

Update `protocol.md` only when an item changes a rule a reader could confuse with the primary analysis. Sensitivities are labeled sensitivities. The acceptance test remains the one-gene ΔC rule, which already failed, plus the external funnel comparison, which already failed. New arms do not get to replace a failed prespecified test.

## What you should be able to say when this is done

Patient identity was corrected, and the shortlist does not depend on which of the four duplicated patients’ aliquots you keep, or it does and you have shown how. The GEO reproducibility gate was checked against an unverified pairing. On the same held-out patients, the funnel’s selected gene does not add the prespecified discrimination, a published score was measured in the same way, and a ridge on the funnel’s own genes is a different and already-reported result. A second null scrambles outcomes without refitting differential expression. The intervals, the event counts, and the testing families are written down as they were computed.
