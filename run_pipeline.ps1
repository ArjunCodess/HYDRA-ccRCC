param(
  [switch]$ForceDownload,
  [switch]$SkipInstall,
  [ValidateRange(1, 5)][int]$NestedWorkers = 1
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

function Invoke-NestedCV {
  param([string]$RscriptPath, [int]$Workers)
  if ($Workers -eq 1) {
    Invoke-RStep -RscriptPath $RscriptPath -Step "analysis/22_nested_cv.R"
    return
  }

  $processes = @()
  try {
    foreach ($worker in 0..($Workers - 1)) {
      $repeatIds = @(1..10 | Where-Object { (($_ - 1) % $Workers) -eq $worker })
      $env:HYDRA_REPEAT_IDS = $repeatIds -join ','
      $stdout = "data/processed/nested_worker_$worker.stdout.log"
      $stderr = "data/processed/nested_worker_$worker.stderr.log"
      $processes += Start-Process -FilePath $RscriptPath -ArgumentList "analysis/22_nested_cv.R" `
        -WorkingDirectory (Get-Location).Path -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr
      Write-Host "Nested worker $worker started for repeats $env:HYDRA_REPEAT_IDS"
    }
  }
  finally {
    Remove-Item Env:HYDRA_REPEAT_IDS -ErrorAction SilentlyContinue
  }

  $failed = @()
  foreach ($process in $processes) {
    $process.WaitForExit()
    if ($process.ExitCode -ne 0) {
      $failed += $process.Id
    }
  }
  if ($failed.Count -gt 0) {
    throw "Nested CV workers $($failed -join ', ') failed. Inspect data/processed/nested_worker_*.stderr.log."
  }
  Invoke-RStep -RscriptPath $RscriptPath -Step "analysis/22_nested_cv.R"
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
  "analysis/00_check_environment.R",
  "analysis/tests/test_identity_and_endpoints.R",
  "analysis/01_download_tcga.R",
  "analysis/02_download_geo_manifest.R",
  "analysis/02_download_geo.R",
  "analysis/00_verify_inputs.R",
  "analysis/03_qc_tcga.R",
  "analysis/04_deg_tcga.R",
  "analysis/04b_paired_deg_tcga.R",
  "analysis/05_inspect_geo_metadata.R",
  "analysis/05_deg_geo.R",
  "analysis/06_reproducibility.R",
  "analysis/07_survival_tcga.R",
  "analysis/07b_apeglm_global_survival_sensitivity.R",
  "analysis/08_enrichment_tcga.R",
  "analysis/10_candidate_table.R",
  "analysis/10b_paired_candidates.R",
  "analysis/10c_compare_prior.R",
  "analysis/11_hardening_outputs.R",
  "analysis/13_external_survival_gse29609.R",
  "analysis/14_external_survival_emtab1980.R",
  "analysis/15_cox_bootstrap_uncertainty.R",
  "analysis/16_cv_clinical_increment.R",
  "analysis/17_hpa_cell_source.R",
  "analysis/19_direct_tumor_purity.R",
  "analysis/20_tracerx_multiregion_transportability.R",
  "analysis/21_checkmate025_treatment_interaction.R",
  "analysis/22_nested_cv.R",
  "analysis/28_null_summary.R",
  "analysis/23_survival_shape.R",
  "analysis/24_funnel_ablations.R",
  "analysis/26_acceptance_report.R",
  "analysis/30_update_readme.R",
  "analysis/25_paper_numbers.R",
  "analysis/09_figures_tcga.R",
  "analysis/27_audit_figures.R",
  "analysis/18_write_manifest.R",
  "analysis/12_validate_outputs.R"
)

foreach ($step in $steps) {
  if ($step -eq "analysis/22_nested_cv.R") {
    Invoke-NestedCV -RscriptPath $rscript -Workers $NestedWorkers
  } else {
    Invoke-RStep -RscriptPath $rscript -Step $step
  }
}

Write-Host ""
Write-Host "Pipeline complete."
