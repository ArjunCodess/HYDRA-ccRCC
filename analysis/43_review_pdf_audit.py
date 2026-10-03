"""Check that the compiled review manuscript contains the traced results."""
import csv
import re
import sys
from pathlib import Path

from pypdf import PdfReader

sys.stdout.reconfigure(encoding="utf-8")
root = Path(__file__).resolve().parents[1]
reader = PdfReader(root / "paper/main.pdf")
pages = [page.extract_text() or "" for page in reader.pages]
text = "\n".join(pages).replace("\u2212", "-")
plain = re.sub(r"\s+", " ", text)
assert "Arjun Vijay Prakash" in plain
assert "Hardening Transcriptomic Prognostic Evidence" in plain
assert "Development draft" not in plain and "??" not in plain
for caption in ("Figure 1:", "Figure 2:", "Figure 3:",
                "Complete priority association table", "AI use"):
    assert caption in plain, f"Missing compiled manuscript content: {caption}"

macro_text = (root / "paper/results_macros.tex").read_text(encoding="utf-8")
macros = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{([^{}]*)\}", macro_text))
for name in ("NestedDelta", "NestedCiLow", "NestedCiHigh", "BenchHydraDelta",
             "BenchHydraCiLow", "BenchHydraCiHigh", "BenchRidgeDelta",
             "BenchRidgeCiLow", "BenchRidgeCiHigh"):
    assert macros[name] in plain, f"Compiled PDF lacks current value for {name}"

with (root / "results/tables/high_confidence_candidate_genes.csv").open(
        encoding="utf-8-sig", newline="") as stream:
    genes = [row["symbol"] for row in csv.DictReader(stream)]
table_pages = "".join(pages[i] for i in range(len(pages))
                      if "Complete priority association table" in pages[i])
compact = re.sub(r"\s+", "", table_pages)
assert genes and all(gene in compact for gene in genes), "Priority table is incomplete"
print(f"PASS: {len(pages)} PDF pages; current estimates, three main figures, "
      f"declarations, and all {len(genes)} priority genes present.")
