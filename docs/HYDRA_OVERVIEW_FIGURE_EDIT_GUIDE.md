# Author-maintained overview

The author supplied the updated `paper/figures/hydra_evidence_overview.svg` on 2026-10-04. It retains the original design, adds the missing 1,117-gene gate, updates HYDRA ΔC to +0.0045 with the benchmark conditional interval −0.0072 to +0.0189, separates the external cohorts and composition results, and labels parallel downstream assessments. It remains the first README image and Figure 1 after the abstract.

Stage 45 now reads the maintained SVG and exports it with headless Chrome/Edge/Chromium. It does not redraw, relabel, or restore older artwork. The viewBox fixes the PDF page dimensions, so a three-times SVG export does not change its paper layout. The producer verifies the SVG bytes remain unchanged and records their SHA-256 hash. Both the pipeline and manuscript build invoke it. Set `HYDRA_OVERVIEW_BROWSER` if the browser executable is elsewhere. Temporary print HTML and the isolated browser profile stay in ignored `paper/qa_render/`.

The author applied the essential corrections. Remaining explanatory details belong in the caption and supplement: the benchmark interval conditions on saved predictions; ClearCode34 is a signed-score adaptation; TRACERx numbers summarize different sampling schemes; zero FDR-supported CheckMate interactions does not establish absence of interaction. The image itself is preserved as supplied.
