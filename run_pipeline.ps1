$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = "C:\Users\pcyou\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
$Rscript = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe"

& $Python (Join-Path $Root "scripts\00_prepare_epithelial_matrix.py")
if ($LASTEXITCODE -ne 0) { throw "Step 00 failed." }

& $Rscript (Join-Path $Root "scripts\01_cluster_epithelial.R")
if ($LASTEXITCODE -ne 0) { throw "Step 01 failed." }

& $Rscript (Join-Path $Root "scripts\02_run_infercnv.R")
if ($LASTEXITCODE -ne 0) { throw "Step 02 failed." }

& $Rscript (Join-Path $Root "scripts\03_build_malignant_states.R")
if ($LASTEXITCODE -ne 0) { throw "Step 03 failed." }

& $Rscript (Join-Path $Root "scripts\04_final_kcn_analysis.R")
if ($LASTEXITCODE -ne 0) { throw "Step 04 failed." }

Write-Host "Core pipeline completed. Final files are in outputs\final."
