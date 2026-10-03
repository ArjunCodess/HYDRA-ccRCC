$ErrorActionPreference = "Stop"

$pdflatex = Get-Command pdflatex -ErrorAction SilentlyContinue
$bibtex = Get-Command bibtex -ErrorAction SilentlyContinue
if (-not $pdflatex -or -not $bibtex) {
  throw "pdflatex and bibtex are required to build paper/main.pdf."
}

$paperDir = Join-Path (Get-Location) "paper"
$validationLog = Join-Path (Get-Location) "results/tables/revision_command_log.csv"
function Invoke-PaperCommand {
  param([string]$Executable, [string[]]$Arguments, [string]$Pass)
  & $Executable @Arguments
  $code = $LASTEXITCODE
  [pscustomobject]@{command = "$(Split-Path $Executable -Leaf) $($Arguments -join ' ') [$Pass]"; exit_status = $code; completed_utc = [DateTime]::UtcNow.ToString('o')} |
    Export-Csv -LiteralPath $validationLog -Append -NoTypeInformation
  if ($code -ne 0) { throw "Paper command failed: $Pass" }
}

Invoke-PaperCommand (Get-Command python).Source @('analysis/45_evidence_overview.py') 'overview export'
Write-Host "Building paper/main.pdf"
Push-Location $paperDir
try {
  Invoke-PaperCommand $pdflatex.Source @('-interaction=nonstopmode', '-halt-on-error', 'main.tex') 'first pass'
  Invoke-PaperCommand $bibtex.Source @('main') 'bibliography'
  Invoke-PaperCommand $pdflatex.Source @('-interaction=nonstopmode', '-halt-on-error', 'main.tex') 'second pass'
  Invoke-PaperCommand $pdflatex.Source @('-interaction=nonstopmode', '-halt-on-error', 'main.tex') 'final pass'
  $finalLog = Get-Content -LiteralPath (Join-Path $paperDir 'main.log') -Raw
  if ($finalLog -match 'undefined (references|citations)|Citation.+undefined|Reference.+undefined|Overfull \\hbox') {
    throw 'Final manuscript still has undefined citations/references or overflowing text.'
  }
}
finally {
  Pop-Location
}

Remove-Item -LiteralPath `
  (Join-Path $paperDir "main.aux"), `
  (Join-Path $paperDir "main.bbl"), `
  (Join-Path $paperDir "main.blg"), `
  (Join-Path $paperDir "main.log"), `
  (Join-Path $paperDir "main.out") `
  -Force -ErrorAction SilentlyContinue

Write-Host "Built paper/main.pdf"
[pscustomobject]@{command = "build_paper.ps1"; exit_status = 0; completed_utc = [DateTime]::UtcNow.ToString('o')} |
  Export-Csv -LiteralPath $validationLog -Append -NoTypeInformation

