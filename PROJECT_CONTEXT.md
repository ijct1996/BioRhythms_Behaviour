# Project Context: BioRhythms_Behaviour

> Single source of truth for scientific context. **Core v1 Scripts 1–3 are frozen.** Extended UR (Kent A–G) lives under top-level `Extended/`.

## Scientific question

How does progressive photoperiod lengthening (12:12 → constant light) alter the balance between circadian and ultradian organisation in mouse locomotor activity, and how do circadian and ultradian rhythm bands co-express across the photoperiod gradient?

## Biological background

Mouse behavioural activity recorded at **10 min** sampling (0.167 h steps) over ~5 d (lengthening cohorts) or ~10 d (NR2B LD/DD). Circadian-scale structure (~22–28 h) is characterised by harmonic subtraction (validation layer). **Primary** wavelet analyses use **RAW** activity; harmonic subtraction confirms whether circadian harmonics are present without driving Core downstream period/co-expression statistics.

## Scope: Core v1 vs Extended UR

| | Core v1 (frozen) | Extended UR (E0–E1) |
|---|------------------|---------------------|
| Location | `Analysis/`, `Functions/`, `Config/` | Top-level `Extended/` |
| Cohorts | C57_LP, NR2B_LP, NR2B_LD_DD | Kent A–G pipeline (same / related sheets) |
| Photoperiods | L12–L24 (2 h steps); L0 DD | + projected-dark under LL (L22 ZT22–24) |
| Scripts | 3 (`run_*`) | 7 wrappers → Kent A–G |
| Confirmatory UR | Residual CWT (AnchorOK, ±5% harmonics) | CarryForward Raw UR (±15% SEL_P360) |

**Parallel claims:** Core Script 3 residual-CWT UR and Extended C CarryForward are **not** interchangeable — see `Analysis/METHODS.md` vs `Extended/METHODS_EXTENDED.md`.

## Hypotheses

### Primary (confirmatory)

Longer photoperiods **weaken circadian organisation** and **strengthen ultradian** locomotor rhythms.

### Secondary (confirmatory)

**Multiscale biological rhythm co-expression** in mouse locomotor activity changes systematically across the photoperiod gradient (CR₂₀₋₂₈ vs UR bands).

## Statistical philosophy

- Block-shuffle surrogates for circadian anchor gating (primary); Lomb–Scargle optional QC column (future).
- Reference `Shared/Research_Knowledge.md` for FDR, effect sizes, exploratory vs confirmatory labelling.
- Default harmonic residual when AnchorOK: **SEL_P360** (Selective Min360).
- Core: wavelet on **RAW**; HSub = validation; Script 3 confirmatory UR on residual CWT.
- Extended: RAW ridges gated by SEL_P360 CarryForward (±15%); LL uses projected L22 dark ZT22–24.

## Environment (locked)

| Item | Value |
|------|-------|
| MATLAB | R2025b |
| Code repo | `C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour` |
| Shared library | `C:\Users\User\Dev\Cursor\Research\Shared` (sibling) |
| Results | Outside repo — `user_paths.m` (gitignored) |
| LL (L24) | 24 h light |
| L0 | 0 h light (constant darkness) — no projected-dark in Core v1 |
| Extended LL | Projected dark = L22 **ZT22–24** |

## Pipeline (Core v1 — frozen)

```
Pipeline_Input xlsx
    → Script 1: run_harmonic_subtraction  (01_HarmonicSubtraction/)
    → Script 2: run_behav_wavelet         (02_Wavelet/ + Handoff/)
    → Script 3: run_across_photoperiod    (03_AcrossPhotoperiod_*/)
```

## Pipeline (Extended UR — E1 wrappers)

```
cd Extended → setup_extended_paths
    → A run_extended_hsub
    → B run_extended_wavelet_ridge
    → C run_extended_period_gate          (CarryForward ±15% SEL_P360)
    → D run_extended_ridge_resync         (+ D.5 companions)
    → E run_extended_across_photoperiod
    → F run_extended_phase_events
    → G run_extended_publication_profiles
```

Docs: `Extended/EXTENDED_CONTEXT.md`, `Extended/METHODS_EXTENDED.md`, `Extended/README.md`.

## Current progress

| Milestone | Status |
|-----------|--------|
| Core v1 folder skeleton | Complete |
| Script 1–3 implementation | Complete — Core frozen for Extended work |
| Full cohort regression (Phase 6) | Not run |
| Legacy v4–v6 archived | Pending Phase 7 |
| Extended E0 tree + Kent_AG copies | Complete |
| Extended E1 thin wrappers + docs | Complete |
| Extended E2 modularisation | Not started |

## Entry points

| Script | Purpose |
|--------|---------|
| `Analysis/run_harmonic_subtraction.m` | Core HSub + QC workbooks + 4 residuals |
| `Analysis/run_behav_wavelet.m` | Core RAW wavelet, jet scalograms, handoff MAT |
| `Analysis/run_across_photoperiod.m` | Core RAW + HSub-validated UR co-expression |
| `Extended/setup_extended_paths.m` | Extended + Core paths |
| `Extended/Analysis/run_extended_*.m` | Kent A–G launchers |

## Publication record

| Field | Value |
|-------|-------|
| Release tag | core-v1.0 (pending) |
| Zenodo DOI | |
| Commit hash | |
