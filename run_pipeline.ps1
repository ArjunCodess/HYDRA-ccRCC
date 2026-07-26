param(
  [switch]$ForceDownload,
  [switch]$SkipInstall
)

$ErrorActionPreference = "Stop"

function Find-Rscript {
  $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
  if ($cmd) {
    return $cmd.Source
  }

  $candidate = Get-ChildItem -Path "C:\Program Files\R" -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending |
    Select-Object -First 1

  if ($candidate) {
    return $candidate.FullName
  }

  throw "Rscript.exe was not found. Install R for Windows or add Rscript to PATH."
}

function Invoke-RStep {
  param(
    [string]$RscriptPath,
    [string]$Step,
    [string[]]$StepArgs = @()
  )

  Write-Host ""
  Write-Host "==> $Step"

  $previousPreference = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  & $RscriptPath $Step @StepArgs
  $exitCode = $LASTEXITCODE
  $ErrorActionPreference = $previousPreference

  if ($exitCode -ne 0) {
    throw "Pipeline step failed: $Step"
  }
}

$rscript = Find-Rscript

if ($ForceDownload) {
  $env:HYDRA_FORCE_DOWNLOAD = "1"
} else {
  Remove-Item Env:HYDRA_FORCE_DOWNLOAD -ErrorAction SilentlyContinue
}

Write-Host "HYDRA-ccRCC reproducible pipeline"
Write-Host "Rscript: $rscript"
& $rscript --version

$steps = @()
if (-not $SkipInstall) {
  $steps += @("analysis/00_install_packages.R")
}
$steps += @(
  "analysis/01_download_tcga.R",
  "analysis/02_download_geo_manifest.R",
  "analysis/02_download_geo.R",
  "analysis/03_qc_tcga.R",
  "analysis/04_deg_tcga.R",
  "analysis/05_inspect_geo_metadata.R",
  "analysis/05_deg_geo.R",
  "analysis/06_reproducibility.R",
  "analysis/07_survival_tcga.R",
  "analysis/08_enrichment_tcga.R",
  "analysis/10_candidate_table.R",
  "analysis/11_hardening_outputs.R",
  "analysis/13_external_survival_gse29609.R",
  "analysis/14_external_survival_emtab1980.R",
  "analysis/15_selection_stability_tcga.R",
  "analysis/16_cv_clinical_increment.R",
  "analysis/17_hpa_cell_source.R",
  "analysis/18_pipeline_bootstrap_tcga.R",
  "analysis/19_direct_tumor_purity.R",
  "analysis/09_figures_tcga.R",
  "analysis/18_write_manifest.R",
  "analysis/12_validate_outputs.R"
)

foreach ($step in $steps) {
  Invoke-RStep -RscriptPath $rscript -Step $step
}

Write-Host ""
Write-Host "Pipeline complete."
