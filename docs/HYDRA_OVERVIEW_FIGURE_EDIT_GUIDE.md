# Evidence overview corrections and export

The overview at `paper/figures/hydra_evidence_overview.svg` is now the first image in README and Figure 1 immediately after the manuscript abstract. The paper includes its matching vector PDF. The submitted tag preserves the original artwork; the revised figure retains its discovery, evidence-funnel, downstream-assessment, and conclusion structure and palette, with simpler typography and fewer secondary diagnostics.

## Corrections completed

- The six selection gates are 8,534 → 3,323 → 1,186 → 1,117 → 538 → 23, including the previously missing covariate-sensitivity gate. The final set is called priority associations.
- GSE40435/GSE53757 replication belongs within selection. External survival, prediction, and composition are explicitly parallel assessments, not further exclusion gates.
- The one-gene procedure reruns TCGA selection inside training folds with GEO evidence fixed. Its primary conditional patient-bootstrap interval is +0.0045 [−0.0082, +0.0174], in 517 patients and 10 × 5 folds. The separate benchmark-bootstrap interval is not substituted here.
- External direction agrees for 21/22 mapped genes in E-MTAB-1980 and 5/21 in GSE29609. E-MTAB has 12 genes passing both FDR support rules; GSE has no same-direction nominal support, and DDC/TCIRG1 reverse with FDR support. These cohorts were previously inspected.
- Marker adjustment retains FDR support for 17/23 genes, while consensus purity retains direction and FDR support for 23/23. Cell-source context does not establish malignant-cell origin or causality.
- The exploratory ridge benchmark has ΔC +0.0257 [ +0.0116, +0.0414 ]; its Brier-difference interval spans zero. It has no external prediction validation. TRACERx, CheckMate, and secondary comparators remain in the supplement.
- The conclusion identifies a priority shortlist and little one-gene prediction gain, without claiming a validated panel or clinical utility.

## Source and regeneration

Run `python analysis/45_evidence_overview.py` from the repository root. This deterministic producer reads `ledger_funnel_counts.csv`, `tcga_kirc_sample_selection_audit.csv`, primary/benchmark nested summaries, external survival summaries and GSE gene results, and both composition sensitivity tables under `results/tables/`. It asserts the serial ordering, cohort reversal identities, candidate denominators, primary/benchmark point-estimate agreement, failed primary prediction criterion, and ridge Brier interval before export. The SVG metadata records SHA-256 hashes of those inputs.

Python 3.14.0 and ReportLab 4.4.7 produce matching SVG/PDF vector geometry and text. ReportLab's invariant PDF mode removes timestamp variability. The producer can use the existing bundled Codex Python if the default Python lacks ReportLab; it does not install dependencies. `run_pipeline.ps1` invokes it after stage 40, and `build_paper.ps1` exports it again before compilation. Stage 41 asserts that this is the first paper figure and README image; stage 43 checks the compiled caption precedes the introduction. Stages 42, 18, and 12 cover its producer, hashes, and required outputs.
