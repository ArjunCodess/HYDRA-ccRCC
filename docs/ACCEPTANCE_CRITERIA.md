# HYDRA-ccRCC Acceptance Criteria

This file defines engineering and evidence gates. It does not supply student submission prose.

## Current status

- The E-MTAB-1980 frozen-candidate test is implemented, has run successfully, and passes automated output validation.
- Thirteen candidates currently pass the revised E-MTAB-1980 same-direction, unadjusted-FDR, and limited-adjustment-FDR rule; proportional-hazards tests are reported as diagnostics rather than exclusion gates.
- Coefficient uncertainty refits the adjusted Cox coefficient for all 27 candidates in 1,000 event-stratified patient bootstraps and passes automated validation.
- Repeated held-out clinical-increment testing runs 20 repeats of five-fold cross-validation and passes automated validation.
- Human Protein Atlas v25.1 cell-source triangulation maps all 27 revised candidates and passes automated validation.
- TRACERx Renal multiregion sensitivity maps all 27 candidates, quantifies regional readout discordance, and completes 1,000 one-region-per-patient repeats in full and 39-patient scenarios.
- CheckMate 025 treatment-interaction analysis maps all 27 candidates in 250 RNA-profiled randomized patients and retains the null multiplicity-controlled result for both survival endpoints.
- Source provenance, input checksums, package versions, random seeds, and output-manifest checksums are recorded.
- The prespecified computational pipeline is complete. Published consensus tumor-purity sensitivity is implemented and passes automated validation.

## Gate 1: Reproducibility

- A clean machine can run the full pipeline from public accession identifiers.
- Downloads are cached, checksummed, and tied to source URLs and access dates.
- Random seeds, package versions, thresholds, and cohort inclusion rules are recorded.
- Generated tables are validated for required columns, row counts, monotonic funnels, and impossible values.

## Gate 2: Frozen external survival validation

- The 27 revised candidates are defined without using E-MTAB-1980 outcomes, but earlier project versions had inspected the cohort, so the current analysis is not described as prospectively frozen or blind.
- E-MTAB-1980 is used only for testing; no result from it changes the candidate definition or TCGA ranking.
- Every candidate is reported, including genes that are missing, non-significant, or directionally contradictory.
- The primary external result is a continuous-expression Cox model. A limited age, T-stage, and metastatic-status adjusted model is secondary because E-MTAB-1980 has only 23 deaths.
- Multiple-testing-adjusted results and raw effect directions are both retained.

## Gate 3: Selection stability and optimism

- [x] Adjusted-coefficient uncertainty is estimated by refitting all 27 candidates in 1,000 patient bootstraps without turning bootstrap significance into another selection gate.
- [x] Predictive value is evaluated out of sample against a clinical-only model.
- [x] Any apparent gain in concordance is reported with uncertainty and without calling the model clinically useful.
- [x] The analysis distinguishes association stability from external validation.

## Gate 4: Composition and biological source

- [x] Published consensus tumor-purity estimates and prespecified immune, endothelial, stromal, and proximal-tubule marker scores are assessed with documented methods.
- [x] Candidate attenuation after composition adjustment is reported as evidence, not treated as a failed analysis.
- [x] Human Protein Atlas v25.1 normal-tissue single-cell data are used to ask which cell compartments express the frozen candidates.
- [x] No RNA-only result is described as a therapeutic target or mechanism.

## Gate 5: Claim discipline

- The project tests the evidence-hardening framework, not a predetermined five-gene success story.
- The abstract and conclusion match the external results, including a null or contradictory result.
- “Candidate prognostic association” remains the strongest gene-level label unless independent outcome validation supports more.
- All citations are verified by the student against the original source.

## Gate 6: IRIS package

- The student independently writes the synopsis, research paper, abstract, poster, citations, and video script.
- A mentor reviews statistics and biological interpretation without taking over the student's work.
- No school, city, or state appears in blinded submission material.
- The source-data age question and the final 2026–27 deadline are confirmed in writing with IRIS.
- The code repository, research notebook, contribution log, and required forms are ready for audit.
