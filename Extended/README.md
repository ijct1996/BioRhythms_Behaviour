# Extended UR — Scripts 4–7 (development)

Extended pipeline beside **frozen Core Scripts 1–3**. Reuses Core `Handoff/` outputs; adds ridge validation, transition resync, LME/FDR, and phase/profile analyses only.

**Plot mode:** `development` (96 DPI PNG) by default. Set `cfg.plotMode = 'publication'` before a later publication pass.

**Primary ultradian bands:** `UR_1_3` and `UR_3_6` (co-primary for transition/resync figures).

## Prerequisites

Core Scripts 1–3 complete for your cohort, with `Handoff/CoreSummary__*.mat` files present.

## Run order (MATLAB R2025b)

```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour\Extended')
setup_extended_paths

% 4 — ridge handoff from Core + CarryForward gate (±15% SEL_P360)
run_extended_script4_ridge_validation   % pick Core Handoff/ folder

% 5 — transition resync, FDR, LL projected, photoperiod gradient
run_extended_script5_transition_resync  % pick ExtendedHandoff/AcrossPhotoperiod_Input/

% 6 — across-photoperiod LME / FDR (local Extended; tables primary in dev)
run_extended_script6_across_lme

% 7 — phase events + profiles (Kent F+G; skipped in dev unless legacy figs enabled)
run_extended_script7_phase_profiles
```

### One-liner paths (non-interactive)

```matlab
handoff = 'D:\Results\C57_LP\Handoff';
extIn   = 'D:\Results\C57_LP\ExtendedHandoff\AcrossPhotoperiod_Input';
run_extended_script4_ridge_validation(handoff)
run_extended_script5_transition_resync(extIn)
run_extended_script6_across_lme(extIn)
run_extended_script7_phase_profiles(extIn)
```

## Script map

| Extended | Replaces Kent | Input | Key outputs |
|----------|---------------|-------|-------------|
| **4** | B ridge + C | Core `Handoff/` | `ExtendedHandoff/AcrossPhotoperiod_Input/WP_*`, `HSubSupported_PeriodMap.mat` |
| **5** | D + D.5 | Script 4 handoff | `Ultradian_RidgePhase_Resync/`, `TransitionEffect_vs_Photoperiod/` |
| **6** | E | Script 4 handoff | LME/FDR workbooks (`AcrossPhotoperiod_LME/`) |
| **7** | F + G | Script 4 handoff | Phase events + profiles (legacy; off in dev by default) |

## What Core already covers (do not rerun)

| Core | Extended does **not** duplicate |
|------|----------------------------------|
| Script 1 | HSub, SEL_P360 residuals |
| Script 2 | RAW scalograms, peaks, band power handoff |
| Script 3 | Co-expression, HSub-validated residual UR, stitched scalograms |

## Development vs publication

```matlab
cfg = extended_defaults();
cfg.plotMode = 'development';   % default — 96 DPI PNG
cfg.plot.generateLegacyFigures = false;  % Script 7 F/G off in dev
extended_script5_run(extIn, cfg);
```

For publication later: `cfg.plotMode = 'publication'; cfg.plot.generateLegacyFigures = true;`

## Deprecated wrappers

`run_extended_hsub`, `run_extended_wavelet_ridge`, `run_extended_period_gate`, etc. redirect to Scripts 4–7 or Core Script 1.
