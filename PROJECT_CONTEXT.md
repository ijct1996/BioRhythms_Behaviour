# Project Context: BioRhythms_Behaviour (Core v1)

> Single source of truth for scientific context. Core cohort only — **no Kent Extended UR (A–G)** in v1.

## Scientific question

How does progressive photoperiod lengthening (12:12 → constant light) alter the balance between circadian and ultradian organisation in mouse locomotor activity, and how do circadian and ultradian rhythm bands co-express across the photoperiod gradient?

## Biological background

Mouse behavioural activity recorded at **10 min** sampling (0.167 h steps) over ~5 d (lengthening cohorts) or ~10 d (NR2B LD/DD). Circadian-scale structure (~22–28 h) is characterised by harmonic subtraction (validation layer). **Primary** wavelet analyses use **RAW** activity; harmonic subtraction confirms whether circadian harmonics are present without driving downstream period/co-expression statistics.

## Scope: Core v1 vs Extended UR

| | Core v1 | Extended UR (later) |
|---|---------|---------------------|
| Cohorts | C57_LP, NR2B_LP, NR2B_LD_DD | Kent A–G sheets |
| Photoperiods | L12–L24 (2 h steps); L0 DD | Projected-dark variants |
| Scripts | 3 (`run_*`) | TBD |

## Hypotheses

### Primary (confirmatory)

Longer photoperiods **weaken circadian organisation** and **strengthen ultradian** locomotor rhythms.

### Secondary (confirmatory)

**Multiscale biological rhythm co-expression** in mouse locomotor activity changes systematically across the photoperiod gradient (CR₂₀₋₂₈ vs UR bands).

## Statistical philosophy

- Block-shuffle surrogates for circadian anchor gating (primary); Lomb–Scargle optional QC column (future).
- Reference `Shared/Research_Knowledge.md` for FDR, effect sizes, exploratory vs confirmatory labelling.
- Default harmonic residual when AnchorOK: **SEL_P360** (Selective Min360).
- Wavelet on **RAW**; HSub = validation only.

## Environment (locked)

| Item | Value |
|------|-------|
| MATLAB | R2025b |
| Code repo | `C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour` |
| Shared library | `C:\Users\User\Dev\Cursor\Research\Shared` (sibling) |
| Results | Outside repo — `user_paths.m` (gitignored) |
| LL (L24) | 24 h light |
| L0 | 0 h light (constant darkness) — no projected-dark in Core v1 |

## Pipeline (Core v1)

```
Pipeline_Input xlsx
    → Script 1: run_harmonic_subtraction  (01_HarmonicSubtraction/)
    → Script 2: run_behav_wavelet         (02_Wavelet/ + Handoff/)
    → Script 3: run_across_photoperiod    (03_AcrossPhotoperiod_*/)
```

## Current progress

| Milestone | Status |
|-----------|--------|
| Core v1 folder skeleton | Complete |
| Script 1 harmonic subtraction | Implemented — MATLAB checkpoint pending |
| Script 2 wavelet + handoff | Implemented — MATLAB checkpoint pending |
| Script 3 across-photoperiod | Implemented — MATLAB checkpoint pending |
| Full cohort regression (Phase 6) | Not run |
| Legacy v4–v6 archived | Pending Phase 7 |

## Entry points

| Script | Purpose |
|--------|---------|
| `Analysis/run_harmonic_subtraction.m` | HSub + QC workbooks + 4 residuals |
| `Analysis/run_behav_wavelet.m` | RAW wavelet, jet scalograms, handoff MAT |
| `Analysis/run_across_photoperiod.m` | RAW + HSub-validated UR co-expression / period comparison |

## Publication record

| Field | Value |
|-------|-------|
| Release tag | core-v1.0 (pending) |
| Zenodo DOI | |
| Commit hash | |
