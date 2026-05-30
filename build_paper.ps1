$ErrorActionPreference = "Stop"

$pdflatex = Get-Command pdflatex -ErrorAction SilentlyContinue
$bibtex = Get-Command bibtex -ErrorAction SilentlyContinue
if (-not $pdflatex -or -not $bibtex) {
  throw "pdflatex and bibtex are required to build paper/main.pdf."
}

$paperDir = Join-Path (Get-Location) "paper"

Write-Host "Building paper/main.pdf"
Push-Location $paperDir
try {
  & $pdflatex.Source -interaction=nonstopmode main.tex
  if ($LASTEXITCODE -ne 0) { throw "pdflatex failed on first pass." }

  & $bibtex.Source main
  if ($LASTEXITCODE -ne 0) { throw "bibtex failed." }

  & $pdflatex.Source -interaction=nonstopmode main.tex
  if ($LASTEXITCODE -ne 0) { throw "pdflatex failed on second pass." }

  & $pdflatex.Source -interaction=nonstopmode main.tex
  if ($LASTEXITCODE -ne 0) { throw "pdflatex failed on final pass." }
}
finally {
  Pop-Location
}

Remove-Item -LiteralPath `
  (Join-Path $paperDir "main.aux"), `
  (Join-Path $paperDir "main.bbl"), `
  (Join-Path $paperDir "main.blg"), `
  (Join-Path $paperDir "main.log") `
  -Force -ErrorAction SilentlyContinue

Write-Host "Built paper/main.pdf"

