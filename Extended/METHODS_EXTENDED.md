# Methods: Extended UR (Kent A–G)

> Pipeline logic for Extended. Core methods remain in root `METHODS.md`. Core Scripts 1–3 are frozen.

## Overview

Extended implements the collaborator **validated-Raw ultradian** pipeline. RAW activity remains the biological signal for rhythm characterisation. Selective harmonic subtraction (**SEL_P360**) is a **validation layer** that decides which Raw UR ridge periods **CarryForward**. Full-ladder residuals are sensitivity / aggressive-subtraction flags only.

## Locked scientific rules

### CarryForward gate (script C)

- Primary HSub mode: **SEL_P360**.
- A Raw UR period candidate carries forward only if a matching Selective-HSub candidate exists within **±15%** relative period tolerance (`PERIOD_TOLERANCE_FRAC = 0.15`).
- Additional ridge coverage / COI QC thresholds apply as coded in Kent C.
- Downstream D–G analyses consume the CarryForward map, not residual-CWT peaks from Core Script 3.

### Projected dark under LL (scripts D, F, G)

- Photoperiod ≥ **24 h** is treated as constant light (LL): no true dark phase.
- Projected reference schedule = preceding **L22** (`LL_PROJECTED_REFERENCE_PHOTOPERIOD_H = 22`).
- Projected dark window for shading / ZT reference = **ZT 22–24** (lights-off of L22 through end of circadian day).
- Projected DL / LD anchors and peri-transition windows are for **reference / aftereffect** analyses only — not literal dark under LL.

### Parallel to Core Script 3 (do not conflate)

| | Core Script 3 | Extended C (+ E/F/G) |
|--|---------------|----------------------|
| Signal for UR claim | HSub **residual** CWT | **Raw** ridge candidates |
| Validation | AnchorOK; exclude harmonics of P₀ (±5%, k≤6) | SEL_P360 period match **±15%** → CarryForward |
| LL dark | Not used in Core v1 | Projected L22 ZT22–24 |

Report both as separate confirmatory strategies if both appear in a manuscript.

## Script roles

| Letter | Role | Entry wrapper |
|--------|------|---------------|
| A | Harmonic subtraction v12 (residuals, recommendations) | `run_extended_hsub` |
| B | Wavelet + ridge / period candidates; `WP_Summary__*.mat` | `run_extended_wavelet_ridge` |
| C | Raw vs SEL_P360 period validation → CarryForward map | `run_extended_period_gate` |
| D | Ultradian ridge phase resync + LL projected aftereffect (FDR) | `run_extended_ridge_resync` |
| D aux / D.5 | Ridge period graphing / plotting companions | noted by D wrapper |
| E | Across-photoperiod validated-Raw UR (FDR, sex-stable LME, PCA fix) | `run_extended_across_photoperiod` |
| F | Validated UR phase–event histograms | `run_extended_phase_events` |
| G | Selected validated UR publication profiles | `run_extended_publication_profiles` |

## Software

| Software | Version |
|----------|---------|
| MATLAB | R2025b |
| Legacy source | Kent Final Versions A–G (copied under `Extended/Legacy/Kent_AG/`) |

## E1 vs E2

- **E1:** wrappers `addpath` Kent_AG and `run('….m')` on monolithic scripts.
- **E2 (in progress):** extract shared logic into `Extended/Functions/{handoff,ridge,stats,plot}/`. Period gate (C) is modular (`extended_period_gate_run`); A/B/D–G may still use Kent wrappers until modularised. See `HANDOFF_SCHEMA_EXTENDED.md`.
