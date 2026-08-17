# Build cohort-level "Supplemental Tables" folder with FigX_* renamed CSV copies.
# Usage: .\build_cohort_supplemental_tables.ps1 -CohortRoot "...\C57_LP"

param(
    [Parameter(Mandatory = $true)]
    [string]$CohortRoot
)

$ErrorActionPreference = 'Stop'
$CohortRoot = (Resolve-Path $CohortRoot).Path
$destRoot = Join-Path $CohortRoot 'Supplemental Tables'
$destCsv = Join-Path $destRoot 'CSV'
New-Item -ItemType Directory -Force -Path $destCsv | Out-Null

$script10Csv = Join-Path $CohortRoot 'Script10_ManuscriptFigures_Publication\Package_Manuscript\Tables\CSV'
if (-not (Test-Path $script10Csv)) {
    $script10Csv = Join-Path $CohortRoot 'Script10_ManuscriptFigures_Publication\Tables\Supplemental\CSV'
}
if (-not (Test-Path $script10Csv)) {
    throw "Script 10 CSV folder not found under $CohortRoot. Run Script 10 first."
}

# SourceFileName -> DestFileName (Fig/panel association)
$map = [ordered]@{
    # Fig 02 — retention + scalogram metadata (Script 14)
    'Script04_Retention_ByBand.csv'              = 'Fig02_E_Retention_ByBand.csv'
    'Script04_Retention_ByPhotoperiod.csv'       = 'Fig02_Supp_Retention_ByPhotoperiod.csv'
    'Script04_Retention_ByPhotoBand.csv'         = 'Fig02_Supp_Retention_ByPhotoBand.csv'
    'Script04_CarryForward_Periods.csv'          = 'Fig02_Supp_CarryForward_Periods.csv'
    'Script04_Matched_Periods_All.csv'           = 'Fig02_Supp_Matched_Periods_All.csv'
    'Script04_QC_Flags.csv'                      = 'Fig02_Supp_QC_Flags.csv'

    # Fig 03 — co-expression
    'AbsolutePower_Summary.csv'                  = 'Fig03_A_AbsolutePower_Summary.csv'
    'CR_UR_Pairs_Summary.csv'                    = 'Fig03_B_CR_UR_Pairs_Summary.csv'
    'LME_Coef_Delta_BH_FDR.csv'                  = 'Fig03_C_LME_Coef_Delta_BH_FDR.csv'
    'LME_Coef_Power_BH_FDR.csv'                  = 'Fig03_C_LME_Coef_Power_BH_FDR.csv'
    'CR_UR_Pairs_PerMouse.csv'                   = 'Fig03_Supp_CR_UR_Pairs_PerMouse.csv'
    'LME_Inference_BH_FDR.csv'                   = 'Fig03_Supp_LME_Inference_BH_FDR.csv'
    'LME_Inference_Raw.csv'                        = 'Fig03_Supp_LME_Inference_Raw.csv'
    'BandConditionSummary_AnalysisUsed.csv'      = 'Fig03_Supp_BandConditionSummary.csv'
    'ValidationMap_Used.csv'                       = 'Fig03_Supp_ValidationMap_Used.csv'
    'Sex_Balance.csv'                              = 'Fig03_Supp_Sex_Balance.csv'
    'Sex_Assignment.csv'                           = 'Fig03_Supp_Sex_Assignment.csv'
    'CR_UR_Pairs_Summary_BySex.csv'              = 'Fig03_Supp_Sex_CR_UR_Pairs_Summary.csv'
    'AbsolutePower_Summary_BySex.csv'            = 'Fig03_Supp_Sex_AbsolutePower_Summary.csv'

    # Figs 4–7 — profile source tables (shared long-format)
    'ActivityComponent_24h.csv'                  = 'Fig04-07_Source_ActivityComponent_24h.csv'
    'PhaseCoherence_24h.csv'                     = 'Fig04-07_Source_PhaseCoherence_24h.csv'
    'ClusterSummary.csv'                           = 'Fig04-07_Source_ClusterSummary.csv'
    'ClusterMembership.csv'                        = 'Fig04-07_Source_ClusterMembership.csv'
    'AllCandidatePeriods.csv'                      = 'Fig04-07_Source_AllCandidatePeriods.csv'
    'PhaseCoherence_DL_LD.csv'                     = 'Fig04-07_Supp_PhaseCoherence_DL_LD.csv'
    'RidgePower_24h.csv'                           = 'Fig04-07_Supp_RidgePower_24h.csv'

    # Figs 5/7 supp — LD Pre/Post R + transition resync audit (Script 5)
    'Resync_PrimaryStats_BH_FDR.csv'             = 'Fig05-07_Supp_Resync_PrimaryStats_BH_FDR.csv'
    'Resync_PrimaryStats.csv'                      = 'Fig05-07_Supp_Resync_PrimaryStats.csv'
    'CandidateDeltaR.csv'                          = 'Fig05-07_Supp_CandidateDeltaR.csv'
    'DeltaR_Summary.csv'                           = 'Fig05-07_Supp_DeltaR_Summary.csv'
    'CandidatePrePost.csv'                         = 'Fig05-07_Supp_CandidatePrePost.csv'
    'PrePostCoherence.csv'                         = 'Fig05-07_Supp_PrePostCoherence.csv'
    'Resync_RealVsPseudoStats_BH_FDR.csv'        = 'Fig05-07_Supp_Resync_RealVsPseudoStats_BH_FDR.csv'
    'BinnedCoherence.csv'                          = 'Fig05-07_Supp_BinnedCoherence.csv'
    'RidgePowerStats_BH_FDR.csv'                 = 'Fig05-07_Supp_RidgePowerStats_BH_FDR.csv'
    'RidgePeriodStats_BH_FDR.csv'                = 'Fig05-07_Supp_RidgePeriodStats_BH_FDR.csv'
    'Resync_ValidatedCandidates_Used.csv'        = 'Fig05-07_Supp_Resync_ValidatedCandidates_Used.csv'
    'LLProj_PrimaryStats_BH_FDR.csv'             = 'Fig05-07_Supp_LLProj_PrimaryStats_BH_FDR.csv'
    'L22_vs_LLProjectedStats_BH_FDR.csv'         = 'Fig05-07_Supp_L22_vs_LLProjectedStats_BH_FDR.csv'
    'LLProj_DeltaR.csv'                            = 'Fig05-07_Supp_LLProj_DeltaR.csv'
    'LL_NoTransitionSummary.csv'                   = 'Fig05-07_Supp_LL_NoTransitionSummary.csv'
    'TransitionEffect_vs_Photoperiod_Summary.csv'  = 'Fig05-07_Supp_TransitionEffect_vs_Photoperiod.csv'

    # Script 9 supplementary figure manifest
    'Script09_Supplementary_FigureManifest.csv'  = 'SuppFig09_FigureManifest.csv'
}

$indexRows = @()
$copied = 0
$missing = @()

foreach ($kv in $map.GetEnumerator()) {
    $srcName = $kv.Key
    $destName = $kv.Value
    $src = Join-Path $script10Csv $srcName
    $dest = Join-Path $destCsv $destName
    if (Test-Path $src) {
        Copy-Item -Path $src -Destination $dest -Force
        $indexRows += [pscustomobject]@{
            DestFile = $destName
            SourceFile = $srcName
            SourcePath = $src
            FigureAssociation = ($destName -replace '_.*$','')
            Status = 'ok'
        }
        $copied++
    } else {
        $missing += $srcName
        $indexRows += [pscustomobject]@{
            DestFile = $destName
            SourceFile = $srcName
            SourcePath = $src
            FigureAssociation = ($destName -replace '_.*$','')
            Status = 'missing_source'
        }
    }
}

# Script 14 scalogram metadata (Fig 02 panels A–D)
$extra = @(
    @{ Src = Join-Path $CohortRoot 'Script14_Scalograms_Publication\FigureManifest_Script14.csv'; Dest = 'Fig02_ABCD_ScalogramManifest.csv' }
    @{ Src = Join-Path $CohortRoot 'Script14_Scalograms_Publication\Settings_Script14.csv'; Dest = 'Fig02_ABCD_ScalogramSettings.csv' }
)
foreach ($e in $extra) {
    if (Test-Path $e.Src) {
        Copy-Item -Path $e.Src -Destination (Join-Path $destCsv $e.Dest) -Force
        $indexRows += [pscustomobject]@{
            DestFile = $e.Dest
            SourceFile = (Split-Path $e.Src -Leaf)
            SourcePath = $e.Src
            FigureAssociation = 'Fig02'
            Status = 'ok'
        }
        $copied++
    }
}

# Script 11 — dominant period supp (cluster context for Figs 4–7)
$s11 = Join-Path $CohortRoot 'Script11_DominantPeriod_Publication\Tables'
if (Test-Path $s11) {
    $s11map = @{
        'DominantPeriod_ByMouse.csv' = 'SuppFig11_DominantPeriod_ByMouse.csv'
        'DominantPeriod_ByCluster.csv' = 'SuppFig11_DominantPeriod_ByCluster.csv'
        'DominantPeriod_PerMouse.csv' = 'SuppFig11_DominantPeriod_PerMouse.csv'
        'DominantPeriod_ClusterSummary.csv' = 'SuppFig11_DominantPeriod_ClusterSummary.csv'
    }
    foreach ($kv in $s11map.GetEnumerator()) {
        $src = Join-Path $s11 $kv.Key
        if (Test-Path $src) {
            Copy-Item -Path $src -Destination (Join-Path $destCsv $kv.Value) -Force
            $indexRows += [pscustomobject]@{
                DestFile = $kv.Value
                SourceFile = $kv.Key
                SourcePath = $src
                FigureAssociation = 'SuppFig11'
                Status = 'ok'
            }
            $copied++
        }
    }
}

# Script 12 — sex-stratified amplitude supp
$s12 = Join-Path $CohortRoot 'Script12_SexStratifiedProfiles_Publication\Tables'
if (Test-Path $s12) {
    Get-ChildItem -Path $s12 -Filter '*.csv' | ForEach-Object {
        $destName = 'SuppFig12_' + $_.Name
        Copy-Item -Path $_.FullName -Destination (Join-Path $destCsv $destName) -Force
        $indexRows += [pscustomobject]@{
            DestFile = $destName
            SourceFile = $_.Name
            SourcePath = $_.FullName
            FigureAssociation = 'SuppFig12'
            Status = 'ok'
        }
        $copied++
    }
}

$indexPath = Join-Path $destRoot 'SUPPLEMENTAL_TABLE_INDEX.csv'
$indexRows | Export-Csv -Path $indexPath -NoTypeInformation -Encoding UTF8

$readme = @"
# Supplemental Tables — $([IO.Path]::GetFileName($CohortRoot))

Consolidated CSV handoff folder. Files are **copies** renamed with figure/panel prefixes.
Original exports remain in Script10/11/12/14 output folders.

## Naming convention

| Prefix | Meaning |
|--------|---------|
| ``Fig02_E_`` | Fig 2 panel E (retention bar) |
| ``Fig02_ABCD_`` | Fig 2 scalogram panels (Script 14) |
| ``Fig02_Supp_`` | Fig 2 / methods QC (Script 4 audit) |
| ``Fig03_A_`` / ``B_`` / ``C_`` | Fig 3 panels A–C source tables |
| ``Fig03_Supp_`` | Fig 3 supplemental / audit tables |
| ``Fig04-07_Source_`` | Long-format profile data for Figs 4–7 |
| ``Fig04-07_Supp_`` | Profile audit tables |
| ``Fig05-07_Supp_`` | LD Pre/Post R + transition resync (Figs 5 & 7 supp bars) |
| ``SuppFig09_`` | Script 9 supplementary figure manifest |
| ``SuppFig11_`` | Script 11 dominant UR periods |
| ``SuppFig12_`` | Script 12 sex-stratified amplitude |

See ``SUPPLEMENTAL_TABLE_INDEX.csv`` for source → dest mapping.

Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')
"@
Set-Content -Path (Join-Path $destRoot 'README.md') -Value $readme -Encoding UTF8

Write-Host "Supplemental Tables folder: $destRoot"
Write-Host "Copied $copied CSV files."
if ($missing.Count -gt 0) {
    Write-Host "Missing $($missing.Count) Script10 CSV sources (skipped): $($missing -join ', ')"
}

# Derived: N mice per cluster × photoperiod (Figs 4–7 panel n)
$pubDir = Join-Path $destRoot 'For_Publication'
$membershipSrc = Join-Path $destCsv 'Fig04-07_Source_ClusterMembership.csv'
if (-not (Test-Path $membershipSrc)) {
    $membershipSrc = Join-Path $pubDir 'Fig04-07_Source_ClusterMembership.csv'
}
$summarySrc = Join-Path $destCsv 'Fig04-07_Source_ClusterSummary.csv'
if (-not (Test-Path $summarySrc)) {
    $summarySrc = Join-Path $pubDir 'Fig04-07_Source_ClusterSummary.csv'
}
if ((Test-Path $membershipSrc) -and (Test-Path $summarySrc)) {
    New-Item -ItemType Directory -Force -Path $pubDir | Out-Null
    $nOut = Join-Path $pubDir 'Fig04-07_Supp_N_ByCluster_Photoperiod.csv'
    $mem = Import-Csv $membershipSrc
    $summary = Import-Csv $summarySrc | Sort-Object BandName, { [int]$_.ClusterRank }
    $pps = @(12, 14, 16, 18, 20, 22, 24)
    $ppLabels = @{12 = 'L12'; 14 = 'L14'; 16 = 'L16'; 18 = 'L18'; 20 = 'L20'; 22 = 'L22'; 24 = 'LL'}
    function Get-SuppBandLabel([string]$bn) {
        switch ($bn) {
            'UR_1_3' { 'UR 1-3' }
            'UR_3_6' { 'UR 3-6' }
            default { ($bn -replace '_', ' ') }
        }
    }
    $nRows = @()
    foreach ($cl in $summary) {
        $sub = $mem | Where-Object { $_.ClusterID -eq $cl.ClusterID }
        $row = [ordered]@{
            Band = Get-SuppBandLabel $cl.BandName
            Cluster = ('C{0:D2}' -f [int]$cl.ClusterRank)
        }
        foreach ($pp in $pps) {
            $n = ($sub | Where-Object { [int]$_.Photoperiod_h -eq $pp } |
                Select-Object -ExpandProperty SignalID -Unique).Count
            $row[$ppLabels[$pp]] = $n
        }
        $nRows += [pscustomobject]$row
    }
    $nRows | Export-Csv -Path $nOut -NoTypeInformation -Encoding UTF8
    Write-Host "Derived N-by-cluster table: $nOut"
}
