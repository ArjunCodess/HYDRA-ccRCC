"""Audit manuscript evidence links, citation keys, generated counts and readability.

Run from the repository root with Python 3. This complements the numerical R
ledger tests; it does not automatically establish causality or scientific validity.
"""
import csv
import hashlib
import re
from pathlib import Path

root = Path.cwd()
main = (root / "paper/main.tex").read_text(encoding="utf-8")
supplement = (root / "paper/supplement.tex").read_text(encoding="utf-8")
all_text = main + "\n" + supplement
bib = (root / "paper/references.bib").read_text(encoding="utf-8")
keys = set(re.findall(r"@\w+\s*\{\s*([^,]+),", bib))
cited = {key.strip() for group in re.findall(r"\\cite\w*\{([^}]+)\}", all_text) for key in group.split(",")}
assert cited <= keys, f"Undefined citation keys: {cited - keys}"
macros = {}
for name in ("results_macros.tex", "review_macros.tex"):
    source = (root / "paper" / name).read_text(encoding="utf-8")
    for key, value in re.findall(r"\\newcommand\{\\(\w+)\}\{([^}]*)\}", source):
        assert key not in macros, f"Duplicate generated macro: {key}"
        macros[key] = value
for key in macros:
    all_text = re.sub(r"\\" + key + r"\b", lambda _: macros[key], all_text)
assert not re.search(r"\\(?:Review\w+|MappedDeg|ReproDeg|MainSurvival|SensitivityPass|StrictGenes|HighGenes)\b", all_text)
for target in re.findall(r"\\(?:input|includegraphics)(?:\[[^]]*\])?\{([^}]+)\}", main + supplement):
    assert (root / "paper" / target).is_file(), f"Missing LaTeX dependency: {target}"
assert main.count(r"\includegraphics") == 3, "The main manuscript must contain exactly three figures."
assert re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", main)[0] == "figures/hydra_evidence_overview.pdf"
assert main.index("figures/hydra_evidence_overview.pdf") < main.index(r"\section{Introduction}")
readme = (root / "README.md").read_text(encoding="utf-8")
assert re.findall(r"!\[[^]]*\]\(([^)]+)\)", readme)[0] == "paper/figures/hydra_evidence_overview.svg"
assert "primary question" in main and "Three secondary questions" in main
assert "validated panel" in main and "Clinical utility is not evaluated" in main
assert "no funding and no conflicts" in main

rows = []
with (root / "paper/claim_evidence.csv").open(encoding="utf-8", newline="") as stream:
    claims = list(csv.DictReader(stream))
assert len({c["claim_id"] for c in claims}) == len(claims)
for claim in claims:
    for source in claim["evidence"].split(";"):
        path = root / source
        assert path.is_file(), f"Missing claim evidence: {source}"
        rows.append({**claim, "source": source, "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "audit_status": "linked; numerical assertions in test_candidate_ledger.R; interpretive claims reviewed manually"})
output = root / "results/tables/claim_evidence_ledger.csv"
with output.open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

# Flag long sentences instead of pretending this heuristic is a biology review.
plain = main[main.index(r"\begin{abstract}"):main.index(r"\bibliographystyle")]
for key, value in sorted(macros.items(), key=lambda item: -len(item[0])):
    plain = re.sub(r"\\" + key + r"\b", lambda _: value, plain)
plain = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}", "", plain, flags=re.S)
plain = re.sub(r"\\(?:cite\w*|ref|label|path)\{[^}]*\}", "", plain)
plain = re.sub(r"\\(?:section\*?|subsection)\{[^}]*\}", "", plain)
plain = re.sub(r"\\\w+\*?(?:\[[^]]*\])?", "", plain)
plain = re.sub(r"[{}$]", "", plain)
sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", plain)
lengths = [len(re.findall(r"\b[\w'-]+\b", text)) for text in sentences if text.strip()]
assert lengths and max(lengths) <= 65, f"Sentence exceeds readability ceiling: {max(lengths)} words"
report = [f"Main prose sentence count: {len(lengths)}", f"Mean words per sentence: {sum(lengths)/len(lengths):.1f}",
          f"Longest sentence: {max(lengths)} words", "Sentences above 40 words require manual inspection:"]
report.extend(text.strip() for text, n in zip([x for x in sentences if x.strip()], lengths) if n > 40)
(root / "results/tables/review_readability.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
print(f"PASS: {len(claims)} linked claims, {len(cited)} citation keys, all LaTeX dependencies, three main figures.")
print(f"Main prose: {len(lengths)} sentences, mean {sum(lengths)/len(lengths):.1f} words, maximum {max(lengths)} words.")
