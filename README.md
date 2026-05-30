# HYDRA-ccRCC

**Hypoxia-Driven Reproducible Analysis of Clear Cell Renal Cell Carcinoma**

**TL;DR:** HYDRA-ccRCC is a reproducible transcriptomics pipeline for testing whether the genes most consistently dysregulated between ccRCC tumor and normal kidney are also the genes that carry prognostic survival signal. It downloads TCGA-KIRC STAR-count RNA-seq data and two GEO validation cohorts, runs DESeq2/limma differential expression, builds a cross-cohort reproducible DEG table, fits adjusted Cox survival models, performs GO enrichment, and generates discordance, evidence-funnel, sensitivity, null-check, and interpretation outputs. The IEEE-style paper is built separately from the analysis pipeline.

HYDRA-ccRCC asks a narrow question: can we move from "this gene is differentially expressed in ccRCC" to "this gene is reproducibly dysregulated and has evidence of survival relevance"? The current answer is useful but disciplined: the first complete run identifies `3,304` reproducible DEGs, `519` strict candidate genes, and `24` high-confidence candidate prognostic associations. The high-confidence set is small enough for focused review, but still requires literature and biological validation before any final gene-level claim.

The research paper lives in [`paper/main.pdf`](paper/main.pdf), with source in [`paper/main.tex`](paper/main.tex) and references in [`paper/references.bib`](paper/references.bib).

## Key Achievements

- **One-command analysis reproducibility:** `.\run_pipeline.ps1` runs package setup, cached data checks/downloads, TCGA QC, TCGA DESeq2, GEO limma validation, reproducibility scoring, Cox survival modeling, GO enrichment, and figure generation.
- **TCGA discovery cohort:** The pipeline downloads TCGA-KIRC STAR-count RNA-seq data and clinical metadata, then analyzes `541` primary tumor and `72` solid tissue normal samples.
- **Independent GEO validation:** The current validation pass includes `GSE40435` with `101` tumor and `101` adjacent non-tumor samples, plus `GSE53757` with `72` tumor and `72` normal samples.
- **Cross-cohort DEG filtering:** TCGA-only tumor-normal differential expression is narrowed into a reproducible evidence table requiring TCGA significance, same-direction GEO support, and nominal validation evidence.
- **Discordance framing:** The central plot compares absolute tumor-normal log2 fold change against absolute adjusted Cox log hazard ratio, directly testing whether tumor-normal magnitude aligns with survival effect.
- **Evidence-funnel framing:** A generated funnel figure shows compression from TCGA DEG to reproducible DEG to survival-associated genes to strict and high-confidence candidates.
- **Interpretation scaffolding:** The pipeline generates a ranked shortlist, survival report, literature-review seed table, cell-type sanity table, threshold sensitivity table, null-overlap check, DEG-versus-prognostic comparison, and generated gene dossiers.
- **Transparent outputs:** Major result tables and figures are written under `results/`, while large downloaded and processed data artifacts are ignored by git.
- **Cautious scientific interpretation:** The current run is pipeline-complete but not publication-final. The broad survival signal is treated as a reason to tighten modeling, not as a finished biomarker list.

## Overview

### What it does

HYDRA-ccRCC prepares a public-data transcriptomics analysis around clear cell renal cell carcinoma. It uses TCGA-KIRC as the discovery cohort, validates tumor-normal expression direction in GEO, and asks whether reproducible tumor-normal genes also show adjusted survival association in TCGA.

The pipeline produces:

- sample summaries
- TCGA PCA figure
- TCGA DESeq2 tumor-normal table
- GEO limma tumor-normal tables
- cross-cohort reproducibility table
- adjusted Cox survival table
- GO Biological Process enrichment table
- candidate evidence table
- ranked high-confidence shortlist
- threshold sensitivity table
- null-overlap check
- DEG-versus-prognostic comparison table
- cell-type sanity table
- literature-review seed table
- generated high-confidence gene dossiers
- volcano and discordance figures
- evidence-funnel figure

### Why it matters

Many beginner cancer transcriptomics projects assume that the largest tumor-normal DEGs are automatically the best findings. HYDRA-ccRCC is built to stress-test that assumption. The project's core scientific idea is:

> Tumor-normal difference is not the same thing as prognostic importance.

That framing keeps the project away from weak "hub gene" overclaims and toward a defensible evidence funnel.

### What is novel here

The novelty is not a new algorithm. The contribution is the analysis discipline:

- require reproducibility across TCGA and GEO
- use continuous-expression Cox models instead of only median-split Kaplan-Meier plots
- compare expression effect size against survival effect size
- interpret findings through a hypoxia-driven ccRCC biology lens
- preserve contradictions instead of cherry-picking clean stories

## Current Results

| Metric | Value |
| --- | ---: |
| TCGA primary tumor samples | 541 |
| TCGA solid tissue normal samples | 72 |
| GSE40435 tumor/normal samples | 101/101 |
| GSE53757 tumor/normal samples | 72/72 |
| TCGA significant genes after symbol mapping | 8,852 |
| Same-direction genes in both GEO cohorts | 9,149 |
| Reproducible DEGs | 3,304 |
| Strict reproducible prognostic candidates | 519 |
| High-confidence candidates | 24 |

The high-confidence count is intentionally conservative: it requires reproducibility, stage/grade-complete Cox significance, proportional-hazards support, meaningful survival effect size, GEO effect support, and same-direction stage/grade sensitivity checks.

The generated shortlist is still not a final biological truth table. It is an auditable worklist for deciding which candidate prognostic associations deserve deeper manual review.

## How To Run

Run the complete cached pipeline from the repository root:

```powershell
.\run_pipeline.ps1
```

Useful options:

```powershell
.\run_pipeline.ps1 -SkipInstall
.\run_pipeline.ps1 -ForceDownload
```

- `-SkipInstall` uses the existing local R package library.
- `-ForceDownload` redownloads TCGA/GEO data instead of using cached processed files.

The latest verified run completed with:

```text
Pipeline complete.
```

If Windows reports that an `.rds` file cannot be replaced, close any open R/RStudio sessions using this project and rerun the command. Interrupted R sessions can keep generated files locked.

## Repository Layout

```text
HYDRA-ccRCC/
  analysis/              ordered R scripts and helper functions
  data/                  generated raw/processed/metadata files, ignored by git
  environment/           sessionInfo.txt plus ignored local R package library
  paper/                 IEEE manuscript source, references, and PDF
  results/               generated tables and figures
  protocol.md            pre-analysis scientific protocol
  run_pipeline.ps1       one-command analysis pipeline entrypoint
  README.md              project overview
```

## Key Files

- [`run_pipeline.ps1`](run_pipeline.ps1): single command that runs the analysis pipeline.
- [`analysis/00_config.R`](analysis/00_config.R): project paths, thresholds, package lists, and dataset settings.
- [`analysis/01_download_tcga.R`](analysis/01_download_tcga.R): TCGA-KIRC STAR-count and clinical metadata download.
- [`analysis/05_deg_geo.R`](analysis/05_deg_geo.R): GEO validation cohort limma analysis.
- [`analysis/06_reproducibility.R`](analysis/06_reproducibility.R): cross-cohort reproducible DEG scoring.
- [`analysis/07_survival_tcga.R`](analysis/07_survival_tcga.R): adjusted Cox survival modeling.
- [`analysis/10_candidate_table.R`](analysis/10_candidate_table.R): final joined evidence table.
- [`analysis/11_hardening_outputs.R`](analysis/11_hardening_outputs.R): ranked shortlist, sensitivity/null checks, cell-type sanity table, literature-review seed table, dossiers, and evidence-funnel figure.
- [`analysis/12_validate_outputs.R`](analysis/12_validate_outputs.R): final output validation step.
- [`protocol.md`](protocol.md): locked research question, thresholds, and limitations.
- [`paper/main.pdf`](paper/main.pdf): first IEEE-style manuscript draft.

## Dataset Rules

Use scripted downloads only. Do not manually edit downloaded metadata in Excel.

- TCGA-KIRC: downloaded from GDC using `TCGAbiolinks` with `workflow.type = "STAR - Counts"`.
- GEO validation: `GSE40435` and `GSE53757` are downloaded as Series Matrix plus supplementary files.
- Raw FASTQ/BAM files are intentionally not downloaded for v1.
- GEO2R result exports are not used as primary data.

## Scientific Guardrails

- Kaplan-Meier plots are secondary visualization only, not the main survival test.
- Cox models use continuous expression and available clinical covariates.
- High-confidence genes are candidate prognostic associations, not candidate biomarkers.
- Adjacent normal tissue is not treated as perfectly healthy tissue.
- Bulk RNA-seq cannot resolve cell-type source.
- Association is not causation.
- No random forest, pan-cancer expansion, or large black-box signature in v1.
- Gene findings should be described as candidate prognostic associations, not therapeutic targets.

## Completion State

The implemented v1 pipeline now includes the stress tests that were previously listed as next priorities:

1. Survival modeling is restricted to reproducible DEGs when the reproducibility table exists.
2. Cox warnings are captured in the survival output table instead of being hidden in console output.
3. The main model uses stage/grade-complete cases, with separate stage-complete and grade-complete sensitivity models.
4. Candidate filtering requires proportional-hazards support, meaningful survival effect size, GEO effect support, and same-direction sensitivity evidence.
5. Candidate tables include pathway-class annotations for hypoxia, angiogenesis, ECM/EMT, metabolism, immune programs, and unclassified genes.
6. The manuscript has been updated to report the compressed candidate set.
7. The pipeline generates evidence-funnel, threshold-sensitivity, null-overlap, DEG-versus-prognostic, literature-seed, cell-type sanity, survival-report, and dossier outputs.

## Remaining Scientific Work

The codebase is pipeline-complete, but the science is not publication-final. The next work is manual candidate-level interpretation: use the generated 24-gene dossiers and literature table to review ccRCC literature, kidney cell-type expression, and independent validation options before making any gene-level claims.

