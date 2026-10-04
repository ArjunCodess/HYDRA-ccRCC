"""Export the author-maintained overview SVG without modifying its artwork."""
import hashlib
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "paper/figures/hydra_evidence_overview.svg"
output = source.with_suffix(".pdf")
original = source.read_bytes()
root = ET.fromstring(original)
_, _, width, height = map(float, root.attrib["viewBox"].split())
assert width > 0 and height > 0
browser = os.environ.get("HYDRA_OVERVIEW_BROWSER") or next((str(p) for p in (
    Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
    Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
) if p.is_file()), None) or shutil.which("chromium") or shutil.which("google-chrome")
if not browser:
    raise SystemExit("Chrome/Edge/Chromium is required to export the maintained SVG; set HYDRA_OVERVIEW_BROWSER.")
qa = ROOT / "paper/qa_render"
qa.mkdir(parents=True, exist_ok=True)
export = qa / "overview_export.pdf"
export.unlink(missing_ok=True)
html = qa / "overview_print.html"
html.write_text(f'''<!doctype html><meta charset="utf-8"><title>HYDRA evidence overview</title>
<style>@page {{ size: {width}px {height}px; margin: 0; }}
html,body {{ margin:0; width:{width}px; height:{height}px; }}
img {{ display:block; width:{width}px; height:{height}px; }}</style>
<img src="{source.as_uri()}" alt="HYDRA evidence overview">''', encoding="utf-8")
subprocess.run([browser, "--headless", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
                f"--user-data-dir={qa / 'overview_browser_profile'}", "--no-pdf-header-footer",
                "--virtual-time-budget=2000", f"--print-to-pdf={export}", html.as_uri()],
               check=True, timeout=55, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
               creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
assert export.is_file() and export.stat().st_size > 1000
assert source.read_bytes() == original, "Export must not modify the author's SVG."
export.replace(output)
print("PASS: author SVG preserved; exported PDF; source SHA-256 " + hashlib.sha256(original).hexdigest())
