# Methods: BioRhythms Core v1

> Pipeline logic (not code). Three scripts, RAW-primary wavelet rule.

## Overview

Pre-split Excel files (one photoperiod per file) are processed per cohort. Harmonic subtraction validates circadian anchor structure. Wavelet scalograms and period peaks use **RAW** activity. Summaries hand off to across-photoperiod co-expression and period comparison plots.

## Pipeline

```
Pipeline_Input/*.xlsx
    ↓  Script 1 — harmonic subtraction (validation + residuals)
01_HarmonicSubtraction/{fileStem}/
    ↓  Script 2 — wavelet on RAW (+ read HSub QC)
02_Wavelet/{fileStem}/  +  Handoff/CoreSummary__{fileStem}.mat
    ↓  Script 3 — cohort-level (Handoff/ only)
03_AcrossPhotoperiod_{cohort}/
```

### Script 1 / 2 batch mode

- Dialog at **start only**: single vs multiple files, folders, mouse columns (Script 1), condition groups (Script 2).
- Then runs **unsupervised** through all selected files (no further prompts).
- Per-file output subfolders named from xlsx stem (`L12_C57_02JUL26`, etc.).
- Run Script 1 then Script 2 on the full cohort before Script 3.

- **Input**: One or more Pipeline_Input xlsx files.
- **Anchor gate**: Regression scan 22–28 h; block-shuffle surrogates (300, 8 h blocks); α = 0.05.
- **Outputs**: Four residuals (Selective/FullLadder × Min360/Min60), QC workbook, **Recommendation** (default **SEL_P360** when AnchorOK), **HSub scalograms** (two-panel Removed|Residual group averages; optional individual at 150 DPI PNG).
- **Entry**: `run_harmonic_subtraction`

## Script 2: Wavelet

- **Primary signal**: RAW activity columns.
- **HSub layer**: Reads Script 1 `HarmonicRemoval_Summary.xlsx` for validation metadata only.
- **CWT**: Morlet (`amor`), 60–1590 min; scalograms **jet** colormap (fixed).
- **Peaks**: Top 5 from mean log₁₀ power spectrum per mouse.
- **Bands**: CR₂₀₋₂₈, UR₁₋₃, UR₃₋₆, UR₆₋₁₂, UR₁₂₋₁₈ for co-expression handoff.
- **Entry**: `run_behav_wavelet`

## Script 3: Across photoperiod

- **Input**: `Handoff/` only (all `CoreSummary__*.mat` for cohort) + Script 1 residual exports via `paths.hsub`.
- **RAW layer (exploratory)**: Co-expression and period comparison from Script 2 handoff (`bandPower`, `periodTable` on RAW).
- **HSub-validated layer (confirmatory)**:
  - **AnchorOK mice only**; SEL_P360 residual from Script 1 `TimeSeries/Residual/`.
  - **UR band power** from residual CWT; **CR₂₀₋₂₈** from RAW CWT.
  - **UR period peaks** (1–18 h) on residual; candidates harmonically related to anchor P₀ (±5%, k=1..6) excluded.
  - Tables: `BandPower_hsub_validated.xlsx`, `PeriodPeaks_hsub_validated.xlsx`, `PeriodPeaks_hsub_validation_log.xlsx`, `Coexpression_ratios_hsub_validated.xlsx`.
  - Figures: `Coexpression_CR_UR_hsub_validated`, `Period_comparison_hsub_validated` (combined scatter); **`Period_Distributions/Period_distribution_hsub_validated_{F|M}`** (box + jitter per sex).
- **Stitched scalograms**: RAW + HSub (visual QC).
- **Runs separately**: C57_LP, NR2B_LP, NR2B_LD_DD.
- **Entry**: `run_across_photoperiod`

## Plot modes

Chosen at batch setup via `prompt_plot_mode()` (development default; last choice remembered in `user_paths.m`).

| Figure type | Development | Publication |
|-------------|-------------|-------------|
| Scalograms (RAW, HSub group avg, stitched) | jet, PNG 96 DPI | jet, JPEG 600 DPI |
| HSub individual (optional QC) | — | 150 DPI PNG (fixed) |
| Other Script 3 line/scatter plots | PNG 96 DPI, default colours | Collaborator palette, TNR, 600 DPI JPEG |

## Software

| Software | Version |
|----------|---------|
| MATLAB | R2025b |
| Shared library | Sibling `Research/Shared` |

## Reproducibility

- Random seed: block-shuffle uses MATLAB default RNG per run — set `rng(seed)` at Script 1 start before publication rerun.
- Pause OneDrive sync during batch runs (Phase 6).
