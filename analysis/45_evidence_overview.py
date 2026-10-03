"""Regenerate the overview SVG and matching vector PDF from audited tables.

Retains the submitted overview's discovery/funnel/evidence/conclusion structure
and palette, with six complete gates and parallel downstream assessments.
No estimates are recalculated; all displayed values come from saved analyses.
"""
import csv
import hashlib
import html
import subprocess
import sys
from pathlib import Path

try:
    import reportlab
    from reportlab.pdfgen import canvas
    from reportlab.lib.colors import HexColor
except ImportError:
    runtime = Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
    if runtime.is_file() and runtime.resolve() != Path(sys.executable).resolve():
        raise SystemExit(subprocess.call([str(runtime), __file__]))
    raise SystemExit("Install ReportLab or use the bundled Codex Python runtime to export the overview.")

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
sources = {}


def table(name):
    path = TABLES / (name + ".csv")
    sources[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


funnel = {r["stage"]: int(r["count"]) for r in table("ledger_funnel_counts")}
primary = table("nested_cv_summary")[0]
benchmark = {r["strategy"]: r for r in table("nested_benchmark_summary")}
em = {r["metric"]: int(r["value"]) for r in table("external_survival_emtab1980_summary")}
gse = {r["metric"]: int(r["value"]) for r in table("external_survival_gse29609_summary")}
external = table("external_survival_gse29609")
reversals = sorted(r["symbol"] for r in external if r["external_present"] == "TRUE"
                   and r["external_same_direction"] == "FALSE" and float(r["external_fdr"]) < .05)
marker = table("candidate_clinical_composition_sensitivity")
purity = table("candidate_direct_tumor_purity_sensitivity")
sample_audit = {r["shortLetterCode"]: r for r in table("tcga_kirc_sample_selection_audit")}
marker_n = sum(float(r["composition_adjusted_fdr"]) < .05 for r in marker)
purity_n = sum(r["same_direction_after_purity"] == "TRUE" and float(r["gene_fdr"]) < .05 for r in purity)
stages = ["discovery_de", "reproducible_de", "adjusted_survival", "covariate_sensitivity", "strict_candidate", "priority_shortlist"]
counts = [funnel[s] for s in stages]
assert all(a >= b > 0 for a, b in zip(counts, counts[1:]))
assert len(marker) == len(purity) == counts[-1]
assert reversals == ["DDC", "TCIRG1"]
assert abs(float(primary["mean_delta_c"]) - float(benchmark["hydra"]["mean_delta_c"])) < 1e-12
assert primary["cv_acceptance"] == "FALSE"

W, H = 1200, 880
out = ROOT / "paper/figures"
out.mkdir(parents=True, exist_ok=True)
pdf = canvas.Canvas(str(out / "hydra_evidence_overview.pdf"), pagesize=(W, H), invariant=1)
pdf.setTitle("HYDRA-ccRCC evidence overview")
pdf.setAuthor("Arjun Vijay Prakash")
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       '<title>HYDRA-ccRCC evidence overview</title>',
       '<desc>Six selection gates followed by parallel external-survival, prediction, and composition assessments.</desc>',
       '<metadata>' + html.escape(str(sources)) + '</metadata>']


def rect(x, y, w, h, fill, stroke="#425466", radius=10):
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{fill}" stroke="{stroke}"/>')
    pdf.setFillColor(HexColor(fill)); pdf.setStrokeColor(HexColor(stroke))
    pdf.roundRect(x, H-y-h, w, h, radius, fill=1, stroke=1)


def text(x, y, value, size=20, bold=False, color="#172b4d", center=False):
    value = str(value)
    font = "Helvetica-Bold" if bold else "Helvetica"
    # Fail rather than silently overflow a panel after a changed result label.
    from reportlab.pdfbase.pdfmetrics import stringWidth
    assert stringWidth(value, font, size) <= W-40
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{"bold" if bold else "normal"}" fill="{color}" text-anchor="{"middle" if center else "start"}">{html.escape(value)}</text>')
    pdf.setFont(font, size); pdf.setFillColor(HexColor(color))
    (pdf.drawCentredString if center else pdf.drawString)(x, H-y, value)


def line(x1, y1, x2, y2):
    svg.append(f'<path d="M{x1},{y1} L{x2},{y2}" stroke="#425466" fill="none"/>')
    pdf.setStrokeColor(HexColor("#425466")); pdf.line(x1, H-y1, x2, H-y2)


def delta(row):
    return f'{float(row["mean_delta_c"]):+.4f} [{float(row["patient_bootstrap_ci_low"]):+.4f}, {float(row["patient_bootstrap_ci_high"]):+.4f}]'


rect(0, 0, W, H, "#ffffff", "#ffffff", 0)
text(600, 43, "HYDRA-ccRCC: hardening prognostic evidence", 29, True, center=True)
rect(20, 65, 1160, 73, "#e7f1fb")
text(40, 95, "Discovery: TCGA-KIRC tumor-normal expression", 22, True)
text(40, 122, f'{sample_audit["TP"]["raw_samples"]} tumor aliquots -> {sample_audit["TP"]["selected_samples"]} patients (one tumor each); {sample_audit["NT"]["selected_samples"]} normal-kidney samples', 21)
rect(20, 155, 1160, 205, "#f0e9fa")
text(40, 185, "Six serial selection gates", 23, True)
labels = [("TCGA", "differential", "expression"), ("GEO expression", "replication", "both cohorts"),
          ("Clinically", "adjusted", "survival"), ("Covariate", "sensitivity", "pass"),
          ("Strict", "effect / FDR", "gate"), ("Priority", "associations", "final gate")]
for i, (count, label) in enumerate(zip(counts, labels)):
    x = 35 + 193*i
    rect(x, 205, 165, 118, "#ffffff")
    text(x+82, 236, f"{count:,}", 29, True, center=True)
    for j, s in enumerate(label): text(x+82, 263+23*j, s, 18, center=True)
    if i < 5:
        text(x+179, 263, ">", 23, center=True)
text(40, 345, "Expression replication: GSE40435 + GSE53757. Counts are genes, not patients.", 19)
text(600, 393, "Downstream assessments are parallel; they do not create new exclusion gates", 21, True, center=True)
for x, fill in [(20, "#e7f1fb"), (415, "#fff0dc"), (810, "#e5f4ed")]:
    rect(x, 412, 370, 315, fill)
text(40, 447, "External survival", 23, True)
text(40, 477, "Previously inspected cohorts", 18)
text(40, 514, f'E-MTAB-1980: {em["same_direction_candidates"]}/{em["platform_present_candidates"]} directions agree', 19, True)
text(40, 541, f'{em["strict_external_support_candidates"]} retain both FDR support rules', 19)
text(40, 580, f'GSE29609: {gse["same_direction_candidates"]}/{gse["platform_present_candidates"]} directions agree', 19, True)
text(40, 607, f'{gse["same_direction_nominal_candidates"]} same-direction nominal support', 19)
text(40, 645, "DDC + TCIRG1 reverse with FDR", 19)
text(40, 691, "Cohort-dependent association", 19, True)
text(435, 447, "Held-out prediction", 23, True)
text(435, 477, f'{primary["n_patients"]} patients; {primary["repeats"]} x {primary["folds"]} folds', 19)
text(435, 506, "TCGA selection rerun in training;", 18)
text(435, 530, "GEO replication evidence fixed", 18)
text(435, 568, "One-gene HYDRA change in C-index", 18, True)
text(435, 595, delta(primary), 19)
text(435, 631, "Conditional patient-bootstrap 95% CI", 17)
text(435, 657, "Little gain beyond clinical variables", 18, True)
text(435, 692, "Project prediction criterion fails", 18)
text(830, 447, "Tissue composition", 23, True)
text(830, 492, f'Marker scores: {marker_n}/{len(marker)} retain FDR', 19, True)
text(830, 521, f'Consensus purity: {purity_n}/{len(purity)} retain', 19, True)
text(830, 548, "direction and FDR support", 19)
text(830, 586, "Bulk proxies give different answers", 18)
text(830, 623, "Cell-source evidence is contextual;", 18)
text(830, 650, "it does not prove tumor-cell origin", 18)
text(830, 691, "No causal interpretation", 19, True)
rect(20, 746, 1160, 114, "#fff8d5")
text(40, 776, f'{counts[-1]} priority associations; no validated clinical panel', 23, True)
text(40, 805, f'Exploratory ridge benchmark change in C-index: {delta(benchmark["ridge_eligible"])}', 20)
assert float(benchmark["ridge_eligible"]["patient_bootstrap_brier_ci_low"]) < 0 < float(benchmark["ridge_eligible"]["patient_bootstrap_brier_ci_high"])
text(40, 835, "Ridge Brier interval spans zero; no external prediction validation or clinical utility evaluation.", 19)
svg.append("</svg>")
(out / "hydra_evidence_overview.svg").write_text("\n".join(svg)+"\n", encoding="utf-8")
pdf.showPage(); pdf.save()
print(f"PASS: overview SVG/PDF; six gates {counts}; table-linked replication, prediction and composition.")
print(f"Export runtime: Python {sys.version.split()[0]}, ReportLab {reportlab.Version}; deterministic vector output.")
