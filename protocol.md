# HYDRA-ccRCC Analysis Protocol

## Background

Clear cell renal cell carcinoma is strongly shaped by VHL/HIF biology, hypoxia signaling, angiogenesis, metabolic rewiring, extracellular matrix remodeling, and tumor microenvironment changes. Many genes differ between tumor and adjacent normal kidney, but tumor-normal difference is not automatically prognostic importance.

## Hypothesis

The strongest tumor-normal transcriptional changes in ccRCC are not necessarily the strongest survival-associated changes. Survival will be associated with a narrower hypoxia-driven adaptation program rather than with every reproducible tumor-normal DEG.

## Datasets

- Discovery: TCGA-KIRC RNA-seq STAR raw counts and clinical data from GDC.
- GEO validation for v1: GSE40435 and GSE53757.
- Secondary datasets are deferred until after the v1 candidate set is manually reviewed.

## Inclusion And Exclusion Criteria

- TCGA expression: include primary tumor and solid tissue normal samples for DEG.
- TCGA survival: include primary tumor samples with usable overall survival time/status.
- GEO validation: include human ccRCC tumor and matched/adjacent normal kidney samples.
- Exclude xenograft/tumorgraft samples from main validation.
- Defer additional validation cohorts until after candidate-level interpretation.

## DEG Method

- TCGA: DESeq2 on raw STAR unstranded counts.
- GEO microarray cohorts: limma with platform-aware preprocessing.
- Paired cohorts: include patient-pair blocking/design where metadata supports it.

## Reproducible DEG Definition

A gene is a reproducible DEG if it satisfies all of:

1. TCGA-KIRC FDR < 0.05.
2. TCGA-KIRC absolute log2 fold change >= 1.
3. Same direction of effect in both v1 GEO validation cohorts.
4. Nominal p < 0.05 in at least one GEO validation cohort.
5. Clean gene identifier mapping across datasets.

## Survival Method

- Outcome: overall survival.
- Model: Cox proportional hazards regression.
- Expression: continuous standardized tumor expression, not median split for primary testing.
- Covariates: age, sex, stage, and grade when usable.
- Main prognostic threshold: stage/grade-complete adjusted Cox FDR < 0.05.
- High-confidence threshold: stage/grade-complete adjusted Cox FDR < 0.01, absolute log hazard ratio at least log(1.5), proportional-hazards p >= 0.05, non-trivial GEO effect support, and same-direction nominal support in both sensitivity models.
- Diagnostics: check proportional hazards assumption with `cox.zph`.
- Sensitivity: fit separate stage-complete and grade-complete models in addition to the main stage/grade-complete model.

## Enrichment Method

- ORA for reproducible/prognostic DEG lists.
- Hallmark class annotation for candidate genes when MSigDB access is available through `msigdbr`.
- Primary interpretation anchor: hypoxia-driven adaptation.
- Secondary categories: angiogenesis, metabolism, ECM remodeling, and immune microenvironment.

## Discordance Figure

- x-axis: absolute DEG magnitude, preferably reproducibility-weighted or meta-analysis log2FC.
- y-axis: absolute adjusted Cox log hazard ratio.
- point size: Cox adjusted significance.
- color: pathway category.
- outline/shape: prognostic threshold pass/fail.

## Known Limitations

- Bulk RNA-seq cannot identify cell-type source of expression.
- Adjacent normal tissue is not perfectly healthy tissue.
- Retrospective cohorts can contain confounding.
- Association does not prove causation.
- No wet-lab validation in v1.
