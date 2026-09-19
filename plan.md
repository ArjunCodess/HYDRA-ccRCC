# HYDRA-ccRCC audit remediation PR

## Summary

Build one PR in six short phases. Each phase ends with regenerated outputs and checks before the next begins. The PR will include corrected code, result tables and figures, `plan.md`, and an updated **development** manuscript and PDF. It will report results that fail the proposed hypotheses. Existing external cohorts remain exploratory because they were previously inspected; this PR will not claim untouched external validation.

## Phases

1. **Lock inputs and fix patient identity.** Record hashes and actual retrieval provenance for cached inputs, pin the R environment, and add assertions for sample joins, GEO pairs, and survival-event coding. Select one TCGA tumor sample per patient throughout the primary pipeline using highest raw-count library depth, breaking ties by barcode; apply the same rule to normals. Add a paired tumor–normal DE sensitivity analysis. Keep the prior results available through Git history and summarize changes from them.

2. **Regenerate the corrected discovery analysis.** Rerun DE, survival, candidate selection, all-gene sensitivity, and downstream analyses from the corrected patient table. Remove fixed assertions that require the old 27-candidate result. Bootstrap reference fits and purity comparisons must use the same patients and expression scale as their comparison fits. Validate all joins, cohort counts, FDR calculations, and candidate coverage before interpreting any changed results.

3. **Measure selection optimism.** Move the complete candidate-selection procedure inside patient-level outer CV, with ten repeats of five folds. Fit preprocessing and expression scaling on training patients; keep every patient’s samples together. Compare clinical-only predictions with the training-selected gene procedure using concordance, three-year Brier score, and calibration. If a training fold selects no gene, use its clinical-only prediction and record the event. Run 200 clinical-only null simulations through the same selection procedure. Retain the existing conditional bootstrap as a separate descriptive analysis, clearly labeled as such.

4. **Test the funnel and biological alternatives.** At matched list sizes, compare the complete rule with survival-only, DE-only, and leave-one-GEO-out selection, plus expression-matched control genes. Use the already inspected external cohorts only for exploratory comparisons and report their contradictory findings, including FDR-significant reversals. Add matched-patient purity and marker-score comparisons, time-varying and nonlinear survival checks, and a TRACERx analysis that holds patient subsets fixed while varying region choice. Record unavailable tests rather than substituting normal-tissue HPA data for tumor single-cell evidence.

5. **Rebuild results and manuscript together.** Regenerate every affected table and figure, then update `README.md`, `protocol.md`, acceptance criteria, and `paper/main.tex` from those outputs. Remove stale counts and unsupported rankings; distinguish observed results from interpretations and proposed work. Preserve null CheckMate interactions, negative overlap findings, composition failures, and external contradictions. Build and inspect `paper/main.pdf`. The paper remains an internal development manuscript, consistent with the repository’s submission guidance.

6. **Verify and open the PR.** Run schema checks, numerical cross-checks, targeted statistical tests, a clean cached-input pipeline run, and the PDF build. Review the diff for generated artifacts and unexplained result changes. Commit the completed work on a `codex/` branch, push it, and open one PR whose description states what changed, which proposed criteria passed or failed, and which validation remains unavailable.

## Acceptance and reporting

The implementation passes its engineering gate only if patient identity and endpoint assertions pass, the pipeline reproduces its committed outputs from pinned inputs, every manuscript number traces to a regenerated artifact, and the PDF builds. Scientific acceptance is separate: the proposed CV increment must be at least 0.01 with a patient-level 95% interval above zero and no worse three-year Brier score; the funnel must outperform its declared survival-only comparator on exploratory external directional replication with a paired 95% interval above zero. Failure of either criterion will be reported as a negative result and will weaken the corresponding manuscript claim; it will not block a truthful PR.

## Assumptions and dependencies

- Use the public data already cached or identified in this repository. New untouched external patients and an independent randomized treatment cohort are future confirmatory work.
- `plan.md` is an untracked user file; include it in the PR without discarding or rewriting its audit findings.
- GitHub CLI is currently unauthenticated. Complete the branch, commits, validation, and PR description first; opening the PR will require a working GitHub authentication method if the current Git credential helper does not provide one.
