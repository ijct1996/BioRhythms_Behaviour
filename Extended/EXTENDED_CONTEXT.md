# Extended UR context (E0–E2)

> Scientific and scope note for the Kent A–G Extended ultradian ridge pipeline. Core v1 remains the frozen primary cohort pipeline.

## Why Extended exists

Core v1 answers photoperiod × circadian/ultradian **co-expression** using RAW wavelet + HSub validation, with Script 3 confirmatory UR from **residual CWT** (AnchorOK, SEL_P360).

Extended UR (Kent A–G) answers a **different** confirmatory path: **validated-Raw ultradian ridges** gated by Selective-HSub period agreement (**CarryForward**), then ridge phase / resync, across-photoperiod FDR analyses, and publication phase/profile figures — including **projected-dark** reference under constant light (LL).

## Scope lock (E0–E2)

| Item | Value |
|------|-------|
| Layout | Top-level `Extended/` only (not mixed into Core `Analysis/`) |
| Code state | E2: C modular; A/B/D–G still E1 Kent wrappers |
| Modularisation | **E2 started** — period gate (C) in `Extended/Functions/handoff/`; A/B/D–G still E1 wrappers |
| Core Scripts 1–3 | **Frozen** — do not change HSub / RAW wavelet / across-photoperiod logic for Extended |
| Default residual label | **SEL_P360** (Selective Min360) |
| CarryForward tolerance | **±15%** period agreement (Kent C `PERIOD_TOLERANCE_FRAC = 0.15`) |
| LL projected dark | Reference **L22**: projected dark **ZT 22–24** when photoperiod ≥ 24 h |

## Parallel claims (do not merge)

| Claim path | What is confirmatory | Gate |
|------------|----------------------|------|
| **Core Script 3** | UR band power / period peaks on **HSub residual CWT**; CR₂₀₋₂₈ from RAW | AnchorOK; harmonic exclusion ±5% of P₀/k |
| **Extended C → E/F/G** | **Raw** UR ridge candidates that **CarryForward** after SEL_P360 match | ±15% period; ridge/COI QC |

These are **parallel**, not sequential. Numbers from one path must not be re-labelled as the other.

## Pipeline (Kent order)

```
A  Harmonic subtraction (v12)     → residuals + QC
B  Behav wavelet (v12)            → RAW + HSub ridges / WP_Summary handoff
C  Raw vs SEL_P360 period gate    → CarryForward map (±15%)
D  Ridge phase resync (+ D.5)     → LL projected L22 anchors
E  Across-photoperiod validated UR → FDR / sex-stable LMEs
F  Phase–event histograms
G  Selected publication profiles
```

## Environment

| Item | Value |
|------|-------|
| MATLAB | R2025b (Kent scripts also note R2025a) |
| Repo | `BioRhythms_Behaviour` |
| Entry | `Extended/setup_extended_paths` then `run_extended_*` |
| Results | Outside repo (dialogs / `user_paths.m`) |

## Progress

| Milestone | Status |
|-----------|--------|
| E0 folder tree + Kent_AG copies | Done |
| E1 thin wrappers + docs + defaults | Done |
| E2 modular refactor of A–G | **C done** (`extended_period_gate_run`); A/B/D–G pending |
| MATLAB end-to-end Extended cohort run | Not started |
