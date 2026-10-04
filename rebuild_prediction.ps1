param([ValidateRange(1, 1)][int]$Workers = 1, [string]$RscriptPath = '')
$ErrorActionPreference = 'Stop'
if (!$RscriptPath) {
    $runtime = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($runtime) { $RscriptPath = $runtime.Source }
    else {
        $runtime = Get-ChildItem -Path 'C:\Program Files\R' -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
        if ($runtime) { $RscriptPath = $runtime.FullName }
    }
}
if (!$RscriptPath -or !(Test-Path -LiteralPath $RscriptPath)) { throw 'Rscript was not found; supply -RscriptPath or add Rscript to PATH.' }
$planSource = @'
source("analysis/00_config.R")
inputs <- c("analysis/31_nested_selection_benchmark.R", "analysis/functions/nested_benchmark.R",
  "analysis/functions/frozen_normalization.R", "analysis/00_config.R", "analysis/functions/patient_samples.R",
  FILES$tcga_se, FILES$tcga_clinical, file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"))
spec <- list(files = unname(tools::md5sum(inputs)),
  screen_n = as.integer(Sys.getenv("HYDRA_SURVIVAL_SCREEN", "2000")), seed = RESAMPLING$seed,
  ridge = "glmnet-5.0-efron-alpha0-lambda.min")
stopifnot(!anyNA(spec$files))
plan <- expand.grid(fold = 1:5, repeat_id = 1:10)
plan <- plan[order(plan$repeat_id, plan$fold), ]
complete <- vapply(seq_len(nrow(plan)), function(i) {
  path <- sprintf("data/processed/nested_benchmark_checkpoints/repeat_%02d_fold_%d.rds", plan$repeat_id[i], plan$fold[i])
  file.exists(path) && identical(readRDS(path)$spec, spec)
}, logical(1))
write.csv(plan[!complete, ], "data/processed/frozen_benchmark_plan.csv", row.names = FALSE)
cat(sum(complete), "matching checkpoints;", sum(!complete), "folds need fitting.\n")
'@
$planSource | Set-Content -LiteralPath data/processed/frozen_benchmark_plan.R -Encoding utf8
& $rscriptPath data/processed/frozen_benchmark_plan.R
if ($LASTEXITCODE -ne 0) { throw 'Benchmark checkpoint planning failed.' }
$plan = @(Import-Csv -LiteralPath data/processed/frozen_benchmark_plan.csv)
$previousSettings = @{}
foreach ($name in @('R_GC_MEM_GROW','R_ENABLE_JIT','HYDRA_REPEAT_IDS','HYDRA_FOLD_IDS')) {
    $previousSettings[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
}
try {
    Remove-Item Env:R_GC_MEM_GROW -ErrorAction SilentlyContinue
    $env:R_ENABLE_JIT = '0'
    foreach ($entry in $plan) {
        $env:HYDRA_REPEAT_IDS = $entry.repeat_id
        $env:HYDRA_FOLD_IDS = (1..([int]$entry.fold)) -join ','
        $label = "repeat_$($entry.repeat_id)_fold_$($entry.fold)"
        Write-Output "Fitting $label in a fresh R process; prior folds replay matching checkpoints."
        $process = Start-Process -FilePath $rscriptPath -ArgumentList 'analysis/31_nested_selection_benchmark.R' `
            -WorkingDirectory (Get-Location).Path -WindowStyle Hidden -PassThru `
            -RedirectStandardOutput "data/processed/frozen_benchmark_$label.stdout.log" `
            -RedirectStandardError "data/processed/frozen_benchmark_$label.stderr.log"
        $process.WaitForExit()
        $process.Refresh()
        $exitStatus = $process.ExitCode
        [pscustomobject]@{command = "R_ENABLE_JIT=0 HYDRA_REPEAT_IDS=$($entry.repeat_id) HYDRA_FOLD_IDS=$env:HYDRA_FOLD_IDS Rscript analysis/31_nested_selection_benchmark.R";
            exit_status = $exitStatus; completed_utc = [DateTime]::UtcNow.ToString('o')} |
            Export-Csv -LiteralPath results/tables/revision_command_log.csv -Append -NoTypeInformation
        if ($exitStatus -ne 0) { throw "Benchmark fold failed: $label." }
    }
    Remove-Item Env:HYDRA_REPEAT_IDS,Env:HYDRA_FOLD_IDS -ErrorAction SilentlyContinue
    & $rscriptPath analysis/31_nested_selection_benchmark.R
    $exitStatus = $LASTEXITCODE
    [pscustomobject]@{command = 'R_ENABLE_JIT=0 Rscript analysis/31_nested_selection_benchmark.R'; exit_status = $exitStatus;
        completed_utc = [DateTime]::UtcNow.ToString('o')} |
        Export-Csv -LiteralPath results/tables/revision_command_log.csv -Append -NoTypeInformation
    if ($exitStatus -ne 0) { throw "Benchmark aggregation failed with exit $exitStatus." }
} finally {
    foreach ($name in $previousSettings.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previousSettings[$name], 'Process')
    }
}
