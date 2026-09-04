# HYDRA-ccRCC Mentor Brief

This is an outreach briefing, not an IRIS/ISEF submission document.

## Project in one paragraph

HYDRA-ccRCC is a reproducible public-data computational oncology project asking whether genes consistently dysregulated between clear cell renal cell carcinoma and adjacent kidney also carry credible prognostic associations. The pipeline analyzes TCGA-KIRC, validates tumor-normal direction in GSE40435 and GSE53757, and stress-tests survival associations across clinical adjustment, composition, external cohorts, multiregion sampling, and a randomized treatment comparison. It narrows 8,852 TCGA-significant genes to 3,304 cross-cohort reproducible DEGs, 536 strict candidates, and 27 high-confidence TCGA-derived associations. GSE29609 is largely contradictory, E-MTAB-1980 supports 13 candidates under the revised rule, and CheckMate 025 supplies no multiplicity-controlled evidence that any candidate modifies nivolumab benefit. These checks keep the result framed as a prioritized set of candidate prognostic associations rather than a validated biomarker panel.

## What has already been built

- Scripted TCGA/GEO acquisition and a one-command R pipeline.
- DESeq2 and limma differential-expression analysis.
- Cross-platform direction-preserving reproducibility filters.
- Continuous-expression Cox models adjusted for age, sex, stage, and grade.
- Proportional-hazards, threshold-sensitivity, null-overlap, and composition screens.
- A complete negative-result-preserving output trail, manuscript draft, and generated figures.
- A verified E-MTAB-1980 test with 101 patients and 23 deaths. Twenty-six candidates mapped, 23 preserved the TCGA direction, and 13 retained same-direction FDR support in both the unadjusted and limited-adjustment models.
- A TRACERx Renal sensitivity that measures regional readout discordance and one-region-per-patient Cox instability for all 27 candidates.
- A CheckMate 025 interaction analysis in 250 RNA-profiled patients. No candidate-by-treatment interaction survived FDR correction for overall or progression-free survival.

## Where expert help matters

### Computational oncology

- Decide whether the renal metabolic-differentiation interpretation is biologically defensible.
- Choose the strongest independent single-cell, spatial, or proteomic triangulation source.
- Separate tumor-cell expression from retained kidney, endothelial, immune, and stromal composition.
- Identify which claims are interesting enough for a conference or journal after IRIS.

### Biostatistics

- Review the candidate-selection procedure for leakage and post-selection bias.
- Design a resampling analysis that repeats the complete selection process.
- Decide how to report external Cox effects with only 23 E-MTAB-1980 deaths.
- Evaluate clinical-model increment without optimistic in-sample concordance claims.
- Review multiple testing, missing covariates, proportional hazards, and calibration.

## The requested mentor commitment

A useful first commitment is one 30–45 minute technical review of the protocol, code, and current outputs, followed by written feedback on the three highest-risk methodological issues. If the fit is good, the student can ask whether the professor or a lab member is willing to advise periodic checkpoints through the IRIS submission. Authorship should be discussed only if later contributions satisfy normal scholarly authorship criteria.

## Materials ready to share

- Git repository and reproducible pipeline.
- `README.md` and `protocol.md`.
- Current manuscript development draft, clearly quarantined from IRIS/ISEF submission.
- Candidate and external-validation tables.
- Acceptance criteria and compliance record.
- A short list of precise statistical and biological questions.
