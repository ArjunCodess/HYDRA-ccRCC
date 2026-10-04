"""Expand the verified dependency map and check maintained artifact coverage.

This records present-file hashes and observed command statuses; it does not
infer successful regeneration from file existence or modification times.
Run from the repository root after the final analysis and manuscript build.
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path.cwd()
DOC = ROOT / "docs/REVIEW_ARTIFACT_SOURCES.md"
PREFIXES = {"T": "results/tables", "F": "results/figures", "P": "paper"}


def expand_braces(value):
    match = re.search(r"\{([^{}]+)\}", value)
    if not match:
        return [value]
    return [item for option in match[1].split(",")
            for item in expand_braces(value[:match.start()] + option + value[match.end():])]


def sha256(path):
    path = ROOT / path
    if not path.is_file():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


producers = {}
descriptions = {}
for line in DOC.read_text(encoding="utf-8").splitlines():
    if not line.startswith("| `"):
        continue
    cells = [cell.strip() for cell in line.split("|")[1:-1]]
    if len(cells) != 3:
        continue
    script = cells[0].strip("`")
    script = script if script.endswith(".ps1") else "analysis/" + script
    descriptions[script] = cells[1]
    directory = None
    for token in re.findall(r"`([^`]+)`", cells[2]):
        if not re.search(r"\.(?:csv|png|pdf|svg|rds|tex|txt|tsv|md|log)$|\.\{", token):
            continue
        if re.match(r"^[TFP]/", token):
            token = PREFIXES[token[0]] + token[1:]
            directory = str(Path(token).parent).replace("\\", "/")
        elif "/" in token:
            directory = str(Path(token).parent).replace("\\", "/")
        elif directory:
            token = directory + "/" + token
        else:
            continue
        for artifact in expand_braces(token):
            paths = ROOT.glob(artifact) if "*" in artifact else [ROOT / artifact]
            for path in paths:
                key = path.relative_to(ROOT).as_posix()
                producers.setdefault(key, []).append(script)

# These stages describe output groups rather than repeating each path in prose.
for name in ("candidate_gene_evidence_table", "strict_candidate_genes",
             "high_confidence_candidate_genes", "candidate_summary"):
    producers.setdefault(f"results/tables/{name}.csv", []).append("analysis/10_candidate_table.R")
producers["results/tables/revision_command_log.csv"] = ["run_pipeline.ps1"]
descriptions["run_pipeline.ps1"] = "Sequential analysis commands and their observed exit statuses."
for name in ("review_artifact_sources.csv", "review_artifact_coverage.txt"):
    producers[f"results/tables/{name}"] = ["analysis/42_review_artifact_manifest.py"]
descriptions["analysis/42_review_artifact_manifest.py"] = "Verified source map, artifact inventory, command log, input and version manifests."

# Every maintained results table/figure must have a producer. Archives are
# historical and deliberately excluded; manuscript source/claim files are inputs.
maintained = {path.relative_to(ROOT).as_posix()
              for directory in ("results/tables", "results/figures")
              for path in (ROOT / directory).rglob("*") if path.is_file()}
maintained |= {"paper/figures/hydra_evidence_overview.pdf", "paper/main.pdf", "paper/results_macros.tex", "paper/review_macros.tex",
               "paper/limitations_macros.tex", "paper/limitations_table.tex",
               "paper/plate_macros.tex",
               "paper/review_candidate_table.tex", "paper/evidence_matrix.tex",
               "paper/benchmark_table.tex", "paper/ridge_spec_table.tex",
               "results/tables/review_artifact_sources.csv", "results/tables/review_artifact_coverage.txt"}
missing = sorted(maintained - set(producers))
assert not missing, "Maintained artifacts lack a producer in the verified map: " + "; ".join(missing)

observed = {}
log = ROOT / "results/tables/revision_command_log.csv"
if log.is_file():
    with log.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            for producer in descriptions:
                if producer in row.get("command", ""):
                    observed[producer] = {key: row.get(key, "") for key in
                                          ("command", "exit_status", "completed_utc")}

rows = []
for artifact in sorted(maintained):
    for producer in sorted(set(producers[artifact])):
        rows.append({"artifact": artifact, "producer": producer,
                     "artifact_sha256": sha256(artifact), "producer_sha256": sha256(producer),
                     "source_description": descriptions.get(producer, "See verified stage dependency map."),
                     "source_map_sha256": sha256(DOC.relative_to(ROOT)),
                     "frozen_input_inventory_sha256": sha256("results/tables/input_manifest.csv"),
                     "version_lock_sha256": sha256("environment/package_versions.csv"),
                     "configuration_sha256": sha256("analysis/00_config.R"),
                     "verification_date": "2026-10-04",
                     "original_input_access_date": "unknown; see source_provenance.csv",
                     "observed_command": observed.get(producer, {}).get("command", ""),
                     "observed_exit_status": observed.get(producer, {}).get("exit_status", ""),
                     "observed_completion_utc": observed.get(producer, {}).get("completed_utc", ""),
                     "execution_caveat": "File hash establishes identity; command status establishes execution. Matching checkpoints may be reused."})

# Avoid a self-referential artifact hash. The final run_manifest supplies the
# manifest file's hash after this writer finishes.
for row in rows:
    if row["artifact"] in ("results/tables/review_artifact_sources.csv", "results/tables/review_artifact_coverage.txt"):
        row["artifact_sha256"] = "self-inventory; final run_manifest records hash"
    elif row["artifact"] == "results/tables/run_manifest.csv":
        row["artifact_sha256"] = "execution manifest rewritten after this inventory; validate its entries with analysis/12_validate_outputs.R"
    elif row["artifact"] == "results/tables/revision_command_log.csv":
        row["artifact_sha256"] = "append-only execution log; changes when subsequent commands complete"
with (ROOT / "results/tables/review_artifact_sources.csv").open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
report = (f"PASS: {len(maintained)} maintained result/manuscript artifacts have verified producers.\n"
          f"Producer edges: {len(rows)}; unmapped maintained artifacts: 0.\n"
          "Historical results/archive files are excluded. Manuscript source and claim files are maintained inputs.\n"
          "Coverage is a dependency check, not proof that all commands have completed.\n"
          "Original input access dates remain unknown; dependency-map verification date is 2026-10-04.\n")
(ROOT / "results/tables/review_artifact_coverage.txt").write_text(report, encoding="utf-8")
print(report, end="")
