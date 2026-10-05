# Revision validation

Status: the 2026-10-04 limitation extension has completed all 50 corrected benchmark fits, both 200-iteration null audits, both plate sensitivity models, derived artifact regeneration, manuscript compilation, and PDF inspection. Final consistency checks passed. Changes are committed in small groups before the authorized plain push. The execution CSV preserves every observed failure, interruption, retry, and success; partial checkpoints never establish completion.

## Preserved submission and scope

The author confirmed submission commit `1f345735013853f3a3c09a088475fc177329369f`, preserved by annotated tag `icbinb-bio-2026-submitted`. The follow-up started from clean commit `4bc0e75892dc5641286eb269dfde1811e36bd0bc` on `icbinb-review-improvements`. The author authorized short lowercase commits after validation, followed by plain `git push`. No pull request was requested.

## Correction and sensitivity evidence

Stage 46 keeps the shortlist fixed while harmonizing cohort adjustment, describing standardized-coefficient heterogeneity, testing eight deterministic missing-covariate scenarios, stratifying by collection site, and screening accession IDs. Stage 49 adds plate-adjusted and mixed-plate tissue contrasts. These analyses address defined threats without replacing the primary funnel or establishing patient independence, missingness mechanisms, causality, or clinical utility.

Stage 47 reproduced predictions from all 50 original primary folds before correcting held-out normalization. Training expression and selections were preserved; Cox prediction models were refitted. This is an audited test-transform correction, not 50 fresh differential-expression fits. Stage 31 requires all 50 benchmark folds to fit under the corrected transform. Stage 48 checks their current source signatures and produces the before/after comparison and exact SVG edit guide.

Stages 50 and 33 checkpoint null iterations and RNG state. Both completed 200 iterations. Stage 51 verifies a validation-only type-coercion repair for the all-missing permutation results. Stage 52 proves that removing an unused input object and invoking garbage collection were the only source changes for 42 already corrected benchmark results; predictions, fold metadata, and RNG were preserved. No pre-correction result was migrated. Original source and numerical outputs are archived under `results/archive/normalization_20261004`; workstation checkpoint backups remain ignored.

Several heavy fits were interrupted when Windows paging exhausted C: storage. The execution log records their actual unsuccessful statuses. Completed checkpoints were retained. The final launcher uses one fresh R process per missing fold, replays saved preceding-fold RNG states, and restores runtime variables afterward. Earlier fits used `R_GC_MEM_GROW=0` under memory pressure. Aggregation subsequently failed once with a locked-binding error and once with an integer/expression error after disabling JIT. The successful final retry unset the garbage-collection override and used `R_ENABLE_JIT=0`, with unchanged scoring code, predictions, and seeds. These observations do not establish the underlying cause. Closing the author's memory-heavy browser and freeing storage allowed fitting to continue. File presence alone does not establish checkpoint validity.

The corrected HYDRA concordance increment is 0.00468 with a conditional primary bootstrap interval of [-0.00786, 0.01744]. Plate adjustment and the mixed-plate contrast preserve direction and FDR support for all 23 genes, but only 19 retain the original expression-magnitude gate. CRYL1, RBM47, KNTC1, and TNFAIP2 fall below that gate in both contrasts. External age/grade/T-category adjustment reduces E-MTAB-1980 same-direction FDR support to one gene; ACADM retains its direction but loses FDR support after T-category adjustment. These are sensitivity findings, not a newly selected panel.

During the ledger commit, Git reported corrupt loose blob `94d2793ec1179f499262d328f9a7b848b2dc9929`. The working ledger produced exactly that Git object ID. The damaged compressed object was preserved in ignored local storage, then reconstructed from the validated working file without changing its bytes or commit history. `git fsck --full --no-dangling` subsequently exited 0. The cause of the corruption is unknown.

## Artifacts, versions, and validation boundaries

`results/tables/revision_command_log.csv` is the authoritative command/status/time record. Stage 41 generates the readable table below. Stage 42 inventories maintained tables, figures, generated manuscript artifacts, producers, inputs, configuration/package hashes, and observed execution status. Stage 18 writes the final output manifest; stage 12 checks its hashes and numerical consistency. Historical archives and ignored workstation-only checkpoint backups are excluded from the published artifact inventory.

R 4.6.1 matches the 258-package lock. Runtime details are recorded in `environment/review_tool_versions.txt` and stage session information; seeds and resampling settings are in the configuration and producer scripts. The run uses frozen cached public inputs, without claiming fresh-download reproduction. Original retrieval dates remain unknown and are distinguished from execution and source-verification dates.

The author-maintained SVG remains first in README and Figure 1 before the introduction, with its bytes preserved during PDF export. The author supplied the corrected prediction labels on 2026-10-05; the source SVG is used directly, with PDF conversion for LaTeX. The 23-page PDF passed its text audit and visual inspection of pages 1, 2, 4, 6, 18, and 23, covering the first pages, all three main figures, the supplementary harmonized table, and the full priority-gene table. No clipping or layout repairs were needed. Conditional prediction intervals omit training and selection uncertainty. Untouched external cohorts, cross-accession patient links, and prospective clinical evaluation remain unavailable from the supplied artifacts.

## Observed commands and exit statuses

Stage 41 regenerates this table from the execution CSV whenever the claim audit runs. Repeated commands are grouped; `results/tables/revision_command_log.csv` preserves every logged invocation and UTC completion time. Exit −1 denotes an intentional process stop. Failed attempts, including labeled layout-gate failures, are followed by the successful retries documented in the execution history below. All required stages and the final manuscript build/audit passed.

| Command | Observed exit statuses |
| --- | --- |
| `HYDRA_REPEAT_IDS=1,2,3,4,5,6,7,8,9,10 Rscript analysis/31_nested_selection_benchmark.R` | 1 |
| `HYDRA_REPEAT_IDS=1,3,5,7,9 Rscript analysis/31_nested_selection_benchmark.R` | -1 |
| `R_ENABLE_JIT=0 HYDRA_NULL_SIMS=0 Rscript analysis/22_nested_cv.R` | 0 |
| `R_ENABLE_JIT=0 HYDRA_NULL_SIMS=0 Rscript analysis/28_null_summary.R` | 0 |
| `R_ENABLE_JIT=0 HYDRA_NULL_SIMS=0 Rscript analysis/33_survival_permutation_null.R` | 0 |
| `R_ENABLE_JIT=0 HYDRA_NULL_SIMS=0 Rscript analysis/50_clinical_null_refits.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/09_figures_tcga.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/12_validate_outputs.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/18_write_manifest.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/23_survival_shape.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/24_funnel_ablations.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/25_paper_numbers.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/26_acceptance_report.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/27_audit_figures.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/30_update_readme.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/32_published_signature_benchmark.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/34_evidence_display.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/37_central_results.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/38_presentation_figures.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/39_candidate_ledger.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/40_review_figures.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/46_limitations_sensitivity.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/48_normalization_report.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/49_plate_de_sensitivity.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/tests/test_candidate_ledger.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/tests/test_limitations_sensitivity.R` | 0 |
| `R_ENABLE_JIT=0 Rscript analysis/tests/test_plate_sensitivity.R` | 0 |
| `R_GC_MEM_GROW unset; R_ENABLE_JIT=0 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW unset; R_ENABLE_JIT=0 Rscript analysis/49_plate_de_sensitivity.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=10 HYDRA_FOLD_IDS=1 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=10 HYDRA_FOLD_IDS=1,2 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=10 HYDRA_FOLD_IDS=1,2,3 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=10 HYDRA_FOLD_IDS=1,2,3,4 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=10 HYDRA_FOLD_IDS=1,2,3,4,5 Rscript analysis/31_nested_selection_benchmark.R` | 0, 1073807364 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=8 HYDRA_FOLD_IDS=1 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=8 HYDRA_FOLD_IDS=1,2 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=8 HYDRA_FOLD_IDS=1,2,3 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=8 HYDRA_FOLD_IDS=1,2,3,4 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=8 HYDRA_FOLD_IDS=1,2,3,4,5 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=9 HYDRA_FOLD_IDS=1 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=9 HYDRA_FOLD_IDS=1,2 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=9 HYDRA_FOLD_IDS=1,2,3 Rscript analysis/31_nested_selection_benchmark.R` | 0, 1 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=9 HYDRA_FOLD_IDS=1,2,3,4 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 HYDRA_REPEAT_IDS=9 HYDRA_FOLD_IDS=1,2,3,4,5 Rscript analysis/31_nested_selection_benchmark.R` | 0 |
| `R_GC_MEM_GROW=0 R_ENABLE_JIT=0 Rscript analysis/31_nested_selection_benchmark.R` | 1 |
| `R_GC_MEM_GROW=0 Rscript analysis/31_nested_selection_benchmark.R` | 1 |
| `R_GC_MEM_GROW=0 Rscript analysis/46_limitations_sensitivity.R` | 0 |
| `R_GC_MEM_GROW=0 Rscript analysis/52_verify_benchmark_memory_cleanup.R` | 0 |
| `R_GC_MEM_GROW=0 Rscript analysis/tests/test_limitations_sensitivity.R` | 0 |
| `Rscript -e invisible(lapply(c("analysis/12_validate_outputs.R", "analysis/31_nested_selection_benchmark.R", "analysis/44_verify_benchmark_cache_equivalence.R"), parse))` | 0 |
| `Rscript analysis/00_check_environment.R` | 0 |
| `Rscript analysis/00_verify_inputs.R` | 0 |
| `Rscript analysis/01_download_tcga.R` | 0 |
| `Rscript analysis/02_download_geo.R` | 0 |
| `Rscript analysis/02_download_geo_manifest.R` | 0 |
| `Rscript analysis/03_qc_tcga.R` | 0 |
| `Rscript analysis/04_deg_tcga.R` | 0 |
| `Rscript analysis/04b_paired_deg_tcga.R` | 0 |
| `Rscript analysis/05_deg_geo.R` | 0 |
| `Rscript analysis/05_inspect_geo_metadata.R` | 0 |
| `Rscript analysis/06_reproducibility.R` | 0 |
| `Rscript analysis/07_survival_tcga.R` | 0 |
| `Rscript analysis/07b_apeglm_global_survival_sensitivity.R` | 0 |
| `Rscript analysis/08_enrichment_tcga.R` | 0 |
| `Rscript analysis/09_figures_tcga.R` | 0 |
| `Rscript analysis/10_candidate_table.R` | 0 |
| `Rscript analysis/10b_paired_candidates.R` | 0 |
| `Rscript analysis/10c_compare_prior.R` | 0 |
| `Rscript analysis/11_hardening_outputs.R` | 0 |
| `Rscript analysis/12_validate_outputs.R` | 0, 1 |
| `Rscript analysis/13_external_survival_gse29609.R` | 0 |
| `Rscript analysis/14_external_survival_emtab1980.R` | 0 |
| `Rscript analysis/15_cox_bootstrap_uncertainty.R` | 0 |
| `Rscript analysis/16_cv_clinical_increment.R` | 0 |
| `Rscript analysis/17_hpa_cell_source.R` | 0 |
| `Rscript analysis/18_write_manifest.R` | 0 |
| `Rscript analysis/19_direct_tumor_purity.R` | 0 |
| `Rscript analysis/20_tracerx_multiregion_transportability.R` | 0, 1 |
| `Rscript analysis/21_checkmate025_treatment_interaction.R` | 0 |
| `Rscript analysis/22_nested_cv.R` | 0 |
| `Rscript analysis/22_nested_cv.R (HYDRA_REPEAT_IDS=9,10; disjoint fold worker)` | -1 |
| `Rscript analysis/23_survival_shape.R` | 0 |
| `Rscript analysis/24_funnel_ablations.R` | 0 |
| `Rscript analysis/25_paper_numbers.R` | 0 |
| `Rscript analysis/26_acceptance_report.R` | 0 |
| `Rscript analysis/27_audit_figures.R` | 0 |
| `Rscript analysis/28_null_summary.R` | 0 |
| `Rscript analysis/30_update_readme.R` | 0 |
| `Rscript analysis/31_nested_selection_benchmark.R` | -1, 0 |
| `Rscript analysis/32_published_signature_benchmark.R` | 0 |
| `Rscript analysis/33_survival_permutation_null.R` | -1, 0, 1 |
| `Rscript analysis/33_survival_permutation_null.R (independent current-input run)` | 0 |
| `Rscript analysis/33_survival_permutation_null.R (verified type-only checkpoint repair)` | 0 |
| `Rscript analysis/34_evidence_display.R` | 0 |
| `Rscript analysis/35_gse53757_pairing.R` | 0 |
| `Rscript analysis/35_gse53757_pairing.R (after regenerated selection)` | 0 |
| `Rscript analysis/36_aliquot_sensitivity.R` | 1 |
| `Rscript analysis/36_aliquot_sensitivity.R (identical-settings retry)` | 1 |
| `Rscript analysis/36_aliquot_sensitivity.R (retry after storage recovery)` | 0 |
| `Rscript analysis/37_central_results.R` | 0, 1 |
| `Rscript analysis/38_presentation_figures.R` | 0, 1 |
| `Rscript analysis/39_candidate_ledger.R` | 0 |
| `Rscript analysis/40_review_figures.R` | 0 |
| `Rscript analysis/44_verify_benchmark_cache_equivalence.R` | 0 |
| `Rscript analysis/44_verify_benchmark_cache_equivalence.R --audit-only` | 0 |
| `Rscript analysis/46_limitations_sensitivity.R` | 0, 1 |
| `Rscript analysis/47_correct_primary_normalization.R` | 0 |
| `Rscript analysis/48_normalization_report.R [initial author artwork run]` | 1 |
| `Rscript analysis/48_normalization_report.R [verified LF helper]` | 0 |
| `Rscript analysis/49_plate_de_sensitivity.R` | 1 |
| `Rscript analysis/50_clinical_null_refits.R` | 0 |
| `Rscript analysis/51_verify_permutation_checkpoint.R` | 0 |
| `Rscript analysis/tests/test_candidate_ledger.R` | 0 |
| `Rscript analysis/tests/test_frozen_normalization.R` | 0 |
| `Rscript analysis/tests/test_identity_and_endpoints.R` | 0 |
| `Rscript analysis/tests/test_limitations_sensitivity.R` | 0, 1 |
| `Rscript analysis/tests/test_nested_benchmark.R` | 0 |
| `Rscript analysis/tests/test_null_checkpointing.R` | 0 |
| `Rscript count corrected-source benchmark signatures (26 of 50)` | 0 |
| `Rscript count corrected-source benchmark signatures (initial shell quoting error)` | 1 |
| `Rscript diagnose permutation comparison types (all 200 selections match after type normalization)` | 0 |
| `Rscript diagnose permutation comparison types (missing locked-library setup)` | 1 |
| `Rscript parse limitation and checkpoint sources` | 0 |
| `Rscript parse new analysis sources` | 0 |
| `benchmark worker 1 intentionally stopped to run conservatively while plate sensitivity completes; saved folds retained` | -1 |
| `benchmark workers 2 and 3 intentionally stopped after C-drive paging pressure; completed folds retained` | -1 |
| `benchmark workers intentionally stopped to repartition four workers; completed folds retained` | -1 |
| `bibtex.exe main [bibliography]` | 0 |
| `build_paper.ps1` | 0 |
| `build_paper.ps1 [final layout gate]` | 1 |
| `bundled Python analysis/43_review_pdf_audit.py` | 0 |
| `bundled-python --version and pypdf version; pdftoppm -v` | 0 |
| `bundled-python [initial PDF text inspection; cp1252 console]` | 1 |
| `bundled-python analysis/43_review_pdf_audit.py` | 0, 1 |
| `correction to preceding attempted worker-1 stop: target was already absent; no termination performed` | 0 |
| `git diff --check` | 0 |
| `git fsck --full --no-dangling` | 0 |
| `pdflatex --version; bibtex --version (existing user cache enabled)` | 0 |
| `pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex [final pass]` | 0, 1 |
| `pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex [first pass]` | 0 |
| `pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex [second pass]` | 0 |
| `pdftoppm -f 1 -l 1 -scale-to 1500 -png paper/main.pdf paper/qa_render/final_inspection` | 0 |
| `pdftoppm -f 2 -l 2 -scale-to 1450 -png paper/main.pdf paper/qa_render/author_overview_final` | 0 |
| `pdftoppm -f 2 -l 2 -scale-to 1600 -png paper/main.pdf paper/qa_render/supplied_20261005` | 0 |
| `pdftoppm -f 21 -l 21 -scale-to 1500 -png paper/main.pdf paper/qa_render/final_inspection` | 0 |
| `pdftoppm -f 3 -l 3 -scale-to 1500 -png paper/main.pdf paper/qa_render/final_inspection` | 0 |
| `pdftoppm -f 4 -l 4 -scale-to 1500 -png paper/main.pdf paper/qa_render/final_inspection` | 0 |
| `pdftoppm -f 5 -l 5 -scale-to 1500 -png paper/main.pdf paper/qa_render/final_inspection` | 0 |
| `python - [final artifact and producer SHA-256 verification]` | 0 |
| `python analysis/41_review_claim_audit.py` | 0 |
| `python analysis/42_review_artifact_manifest.py` | 0, 1 |
| `python analysis/45_evidence_overview.py` | 0, 1 |
| `python analysis/45_evidence_overview.py [deterministic regeneration SHA-256 check]` | 0 |
| `python inline artifact/producer/source-map/input-inventory SHA-256 verification` | 0 |
| `python.exe analysis/45_evidence_overview.py [overview export]` | 0 |
| `rebuild_prediction.ps1 -Workers 1 (R_GC_MEM_GROW=0 worker stopped at 35 folds to switch to per-fold processes)` | 1 |
| `rebuild_prediction.ps1 -Workers 1 (default-GC worker stopped to restart with R_GC_MEM_GROW=0)` | 1 |
| `rebuild_prediction.ps1 -Workers 1 (intentionally stopped for C-drive storage)` | 1 |
| `rebuild_prediction.ps1 -Workers 1 (resumed; intentionally stopped at 30 folds for C-drive storage)` | 1 |
| `run_pipeline.ps1 (SkipInstall=True; ForceDownload=False; NestedWorkers=1; StartAt=analysis/18_write_manifest.R)` | 0 |
| `run_pipeline.ps1 (SkipInstall=True; ForceDownload=False; NestedWorkers=1; StartAt=analysis/34_evidence_display.R)` | 0 |
| `run_pipeline.ps1 (SkipInstall=True; ForceDownload=False; NestedWorkers=1; StartAt=analysis/46_limitations_sensitivity.R)` | 0 |

## Author confirmation record

On 2026-10-01 the user supplied the author name **Arjun Vijay Prakash** and reported no funding and no conflicts of interest. These facts are used in the manuscript. The user subsequently confirmed conceptualization, data curation, software, formal analysis, visualization, and writing as the author contribution roles.

OpenAI Codex assisted with the code, evidence audit, figures, and manuscript revision in this session. This records observed assistance, not certification that the author has independently verified the final draft. The study uses public de-identified datasets; no new recruitment, intervention, or institutional ethics approval is claimed.

## Earlier revision execution history

The initial cached-input pipeline regenerated the primary differential-expression, replication, adjusted-survival, and candidate outputs successfully. Their serial counts remained 8,534 → 3,323 → 1,186 → 1,117 → 538 → 23. The GSE53757 pairing sensitivity was repeated after current candidate generation and retained 23/23 genes in the unpaired expression gate.

Stage 36 subsequently stopped with `default method not implemented for type 'expression'`, at `is.finite(partial)` during the summed-counts rule. Its command exit status was 1. An identical-settings retry stopped with `cannot allocate vector of size 147.4 Mb`, also exit 1; C: then had approximately 21 MB free. After the user requested a storage recheck, C: had approximately 19.4 GB free and the stage was restarted. These failures are execution failures, not scientific null results.

The stage-36 repair retries only that exact expression/partial error, at most three times with unchanged DESeq2 settings. Other errors still stop the stage. It neither disables Cook outlier replacement nor changes the design, thresholds, sample rules, or resampling counts. The summary records `deseq_attempts`. A `StartAt` option allows the pipeline to resume after completed earlier stages without refitting them unnecessarily; successful prior-stage output is required before resumption.

After storage recovered, stage 36 completed with exit 0 at 2026-10-01 16:44 UTC. All three rules required one DESeq2 attempt and reproduced priority counts of 19, 20, and 19. The pipeline then resumed at stage 11. The observed storage/memory failures do not justify attributing the earlier expression-type error to a particular internal cause.

The resumed run passed stages 11, 13–17, and 19, then stage 20 stopped with a dplyr `cannot change value of locked binding for 'k'` error while recording direction agreement. That diagnostic-column assignment was replaced with direct vector assignment; sampling and model calculations were preserved. Resuming at stage 20 completed TRACERx and CheckMate regeneration successfully, then entered the nested-fold analysis. This is an execution repair, not a new sensitivity analysis.

The pre-existing 50 primary nested-fold checkpoints do not match all current source signatures. Input-data hashes match. The nested-script difference is explained by CRLF versus LF, and the patient-helper history shows an added aliquot helper, but the recorded configuration hash cannot be recovered from Git history, even after checking line-ending variants. Therefore those checkpoints are not treated as verified current-source reuse: the normal pipeline must refit them. This preserves the existing signature check and may increase execution time substantially. Benchmark checkpoints undergo their own signature checks.

The separate benchmark-cache audit recovered its historical source hashes. Parsed source differs only in a figure output basename, an added inner-fold validation guard, and an unused aliquot helper; input, configuration, seeds, and model specifications match. Stage 44 checks those exact differences and also requires freshly refitted primary-fold predictions, selected genes, patient IDs, and outcomes to agree with each benchmark checkpoint. Its first migration passed with exit 0 for 16 of 50 folds. Originals are retained in `data/processed/benchmark_cache_originals/`, with original and migrated hashes and source references recorded. All model contents remain identical. This is verified checkpoint reuse, not a claim that those benchmark models were newly refitted. Unverified folds retain their old signature and must be refitted unless a later audit establishes the same equivalence.

An auxiliary worker for repeats 9–10 was attempted using the existing `HYDRA_REPEAT_IDS` option. Although available physical memory appeared sufficient, Windows expanded its paging file and free C: storage fell from approximately 10.3 GB to 2.1 GB. That auxiliary worker was intentionally stopped before finishing its first fold (exit −1). It contributes no completed checkpoint or result. The original sequential run continued, retaining its completed folds. This operational attempt does not change folds, seeds, models, or the required 50-fold result.

The fresh primary fold at repeat 4/fold 1 used the recorded no-replacement fallback and selected a different gene from the old benchmark. Repeat 2/fold 5 also differed in replacement policy, even though its saved clinical/HYDRA predictions agreed. Those checkpoints fail the strengthened reuse check. Stage 31 now replays each primary fold's recorded replacement setting for every comparator and retries that same setting on the exact runtime error; it cannot switch settings independently. Stage 44 normalizes only this documented policy replay when comparing source, and requires matching recorded settings as well as fresh primary prediction agreement. Previously migrated metadata can be refreshed only after unchanged model contents and the original backup specification are verified. This repair addresses fair comparator preprocessing; it does not add a selection method or erase the primary run's runtime-dependent fallback limitation.

All 50 fresh primary-fold checkpoints passed their current signature checks. The regenerated primary summary records mean concordance increment 0.0044663, conditional 95% interval −0.0081638 to 0.0174153, mean Brier difference −0.0009153, two replacement fallbacks, and a failed project criterion. The cache audit verified 45 benchmark checkpoints for reuse; five preprocessing-policy mismatches require refitting. These values supersede the submitted primary estimates rather than being forced to match them. Stage 22 still has its 200 clinical-null simulations to finish before its complete command can return success.

After DE fitting ended, the independent 200-permutation check was started alongside the clinical-null simulations using the five verified fresh first-repeat caches. Both processes use approximately 2–3 GB of private memory, substantially less than the earlier auxiliary DE fit. Stage 33 now saves a complete result checkpoint only after all iterations and result checks pass. Reuse requires matching source/helper/configuration hashes, expression and clinical inputs, all five cache hashes, R/package versions, and resampling count. This is a resampling checkpoint, not newly refitted models at every invocation. A duplicate pipeline invocation must not be allowed to repeat that calculation while the independent producer is still active.

Stage 22 and the clinical-null summary completed with exit 0 at 2026-10-01 20:47 UTC; the selector chose a gene in 3/200 clinical-null simulations. When the pipeline entered stage 33, its duplicate process (PID 43100, confirmed as the child of the task's pipeline PowerShell process and running that exact script) was intentionally stopped. The original independent producer continued. The duplicate records exit −1 and its pipeline wrapper exits 1; neither is a failed scientific null result. Resumption must use the completed, matching-signature permutation checkpoint rather than reduce the 200 iterations or repeat the active calculation.

The independent permutation run completed all 200 iterations with exit 0 at 2026-10-01 21:12 UTC and selected no gene. Its complete checkpoint includes source/input/cache/version signatures and the completion time. After its process exited, C: had approximately 5.9 GB free, physical memory had approximately 9.9 GB available, and commit headroom was approximately 10.4 GB. The pipeline resumed from stage 33 to verify/reuse that result, run the benchmark cache audit, refit the five rejected benchmark folds sequentially, and finish downstream generation and validation.

The resumed stage 33 verified and reused its complete current-input checkpoint with exit 0; stage 44 also passed. Stage 31 refitted repeat 2/fold 2 and saved a current checkpoint matching the primary selected gene. During repeat 2/fold 4, C: fell from approximately 5.9 GB to approximately 70 MB free. The task-owned benchmark process was intentionally stopped before exhausting the drive (stage exit −1, wrapper exit 1); the first completed refit is retained. Four rejected folds remain. More working disk space was requested from the author. This is an observed resource block, not permission to reduce folds, change models, accept stale checkpoints, or report final validation as complete.

`results/tables/revision_command_log.csv` supplies command statuses and UTC completion records. Final completion and artifact validation remain pending until the remaining stages and final manuscript checks finish. Source existence, file timestamps, and earlier successful checks are not substitutes for that completion record.

On the author's 2026-10-03 continuation request, C: had 28.45 GB free. The pipeline resumed sequentially at stage 31 from the 46 current benchmark checkpoints, preserving all earlier completed results. No folds, seeds, thresholds, or model settings were reduced to accommodate storage. Final completion still depends on the remaining refits and downstream checks.

The four resumed benchmark refits completed on 2026-10-03. Repeat 5/fold 1 required an unchanged-settings retry after the known trimmed-mean runtime error and then succeeded. The final audit-only stage-44 command returned 0 and verified all 50 current checkpoints: 45 verified-equivalent reused models and five fresh refits in total. Every current cache matched the fresh primary predictions, outcomes, selected gene, source specification, and recorded replacement setting. Its migration-eligible count of zero means all files already have current specifications, not that any current file failed verification. Stage 31 is still completing its uncertainty summaries before downstream regeneration.

Stage 31 and the published-signature stage 32 completed with exit 0 on 2026-10-03. Stage 37 then stopped with exit 1 because it hard-coded the submitted ridge concordance increment. The current regenerated increment is 0.0256900994, rather than 0.0256882591. That assertion was replaced by recomputing clinical and ridge concordance directly from held-out patient predictions in each repeat, checking all ten repeat deltas, and checking their mean against the summary. Models and result estimates were not changed. The audit therefore verifies internal numerical consistency without forcing a refreshed analysis to reproduce a historical point estimate.

Stage 37 then passed its prediction-based check. Stage 38 stopped with exit 1 because a second legacy assertion expected HYDRA rounded to 0.0039. Its six fixed prediction targets were replaced by strict agreement with the current ten-repeat benchmark and ClearCode34 metric tables. The final prediction validator supplies the direct patient-level agreement checks; presentation no longer enforces historical rounded values.

All downstream writers, ledger tests, claim checks, and artifact-coverage checks completed successfully. Stage 12 then stopped with exit 1 because the newly added patient-prediction checks used dplyr joins without attaching the package. The validator now explicitly loads the already locked and installed dplyr package. No model or result was modified. The manifest is refreshed before rerunning validation because the validator source hash changed.

The final cached-input output validator passed on 2026-10-03. The first refreshed manuscript build completed every native TeX/BibTeX pass but its strict wrapper rejected a 9.29-point overfull abstract line after abbreviation expansion. The cohort sentence was split into two shorter sentences without changing the datasets or interpretation, and the manuscript was rebuilt.

Final completion on 2026-10-03: all required cached-input stages passed, including the direct prediction-agreement and output-hash validator. Splitting the abstract cohort sentence reduced but did not eliminate overflow; shortening the question resolved it. The rebuilt PDF has 21 pages and no undefined citations/references or overfull boxes. Rendered inspection covered pages 1, 3, 4, 5, and 21, including the three main figures and all priority genes. The first ad hoc PDF text probe failed only when printing a Unicode minus through a cp1252 console, after its content assertions passed; the saved UTF-8-aware stage-43 audit returned 0. The funnel was regenerated after inspection to label external results as unadjusted, show marker retention separately, and use current four-decimal prediction estimates. Final inventory/manifest validation follows the completed documentation and PDF.

## Import and manuscript inspection notes

All 37 purity-import warnings concern text `NaN` in the unused IHC column F. The used CPE column produced no such warning. This check does not validate purity biology. The rendered review checked cohort-specific effect scales, the unfavorable ACADM cohort, DDC's separate matched populations, conditional intervals, and the signed ClearCode34 adaptation. Scientific limitations and full methods are in the manuscript and supplement.

## Documentation cleanup

On 2026-10-03, nine obsolete or redundant documents were removed. Unique execution and author-confirmation evidence was consolidated here; claim and README links were updated. Derived tables, figures, ledger, and manuscript were rebuilt from the verified current analysis outputs. Model inputs and specifications did not change, so expensive fitted analyses were preserved rather than refitted for documentation changes.

## Overview placement follow-up (2026-10-03)

The SVG/PDF overview was rebuilt from recorded tables and visually inspected at standalone size and in the compiled paper. It occupies page 2 after the abstract, before the introduction on page 3. README uses the same SVG as its first image. Repeated export preserved byte-identical SVG and PDF hashes. The compiled 22-page manuscript retains three main figures, all 23 priority genes, and current numerical macros. No statistical analysis was refitted for this display change. Stages 41/43 check source and compiled placement; stages 42/18/12 refresh producer coverage and output hashes. Commands and actual exit statuses are appended to the execution log.

## Completion audit follow-up (2026-10-04)

The original task and subsequent figure, document-cleanup, build, and commit requests were checked against maintained artifacts and executable audits. All requested scientific and manuscript changes are implemented. The readable command table was stale after overview inclusion; stage 41 now regenerates it from the execution CSV before hashing claim sources. The requirements audit is recorded in `ICBINB_REVIEW_RESPONSE_MATRIX.md` rather than a new redundant document. Current environment/input locks, patient identity/endpoints, benchmark helpers and cache equivalence, the candidate ledger, claim sources, compiled PDF, producer coverage, and manifest hashes are checked again. No new model fit is justified by this recordkeeping correction. Scientific limitations remain unchanged and explicit in the manuscript.

## Original artwork restoration (2026-10-04)

The author requested the original image rather than the redesign. Both original submitted SVG/PDF files are now restored unchanged; the restoration producer checks exact Git-blob bytes and SHA-256 hashes. First placement in README and the paper is retained. The caption identifies older labels and directs readers to current results. Earlier regeneration statements describe the subsequently reverted redesign. The manuscript is rebuilt and source/PDF/producer/manifest audits are rerun before committing.

## Author-updated overview (2026-10-04)

The author supplied a corrected SVG retaining the original design. Stage 45 now exports that maintained input to PDF with a headless browser, verifies unchanged SVG bytes, and never restores an older version. The updated image stays first in README and after the abstract in the paper. Historical restoration/redesign entries above describe superseded versions. Current numerical results and manuscript audits remain authoritative. The author subsequently authorized building and committing all changes in short lowercase commits.

The final author-image run regenerated derived outputs from stage 34 through stage 27 using the existing source analyses. Stage 42 encountered a transient invalid-argument file-write error and passed on retry. The earlier original-artwork build also had a transient disk-space failure before its successful retry. Failures remain recorded rather than being removed. The author SVG SHA-256 is `eab7a6a7cde3282d22e1e7c490ebb47f6f48b88e366d769f4366b648090358df`; browser export checks it unchanged. Final compilation, claim/PDF checks, artifact coverage, and output hashes are refreshed after regeneration.

## Current author overview, 2026-10-05

The author supplied a new SVG with corrected benchmark labels, the HYDRA interval, unadjusted external-direction wording, and the joint-ridge conclusion. This supplied artwork is used directly as the first README image and converted to PDF for Figure 1. No diagram is generated or redesigned. Earlier artwork warnings describe superseded versions. The caption and README no longer label current values as pre-correction.

The new-artwork build passed, and page 2 was rendered and visually checked with no clipping. Claim and PDF audits passed. The first stage-48 attempt rejected the local normalization helper because Git had converted LF to CRLF; its LF bytes exactly matched the saved MD5. Restoring those line endings made the signature check pass without changing code, checkpoint metadata, fitted models, or results.
