# Extended UR — Scripts 4–12

Extended pipeline beside **frozen Core Scripts 1–3**. Reuses Core `Handoff/` outputs; adds ridge validation, transition resync, LME/FDR, phase/profile analyses, and manuscript figure assembly.

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

% 8 — publication composites (earlier layout; optional once Script 10 is locked)
run_extended_script8_publication_figures   % pick cohort root (e.g. C57_LP)

% 9 — supplementary (sex + cluster grids)
run_extended_script9_supplementary_figures

% 10 — locked manuscript figures (Figs 1–7) + Package_Manuscript handoff zip folder
run_extended_script10_manuscript_figures   % pick cohort root (e.g. C57_LP)

% 11 — dominant UR mean+median periods from Script 7 membership (no WP_TS)
run_extended_script11_dominant_periods     % pick cohort root (e.g. C57_LP)

% 12 — Female|Male activity + post-LD first-peak amplitude (F vs M BH-FDR; after Script 7)
run_extended_script12_sex_stratified_profiles  % pick cohort root (e.g. C57_LP)
```

### One-liner paths (non-interactive)

```matlab
handoff = 'D:\Results\C57_LP\Handoff';
extIn   = 'D:\Results\C57_LP\ExtendedHandoff\AcrossPhotoperiod_Input';
cohort  = 'D:\Results\C57_LP';
run_extended_script4_ridge_validation(handoff)
run_extended_script5_transition_resync(extIn)
run_extended_script6_across_lme(extIn)
run_extended_script7_phase_profiles(extIn)

cfg = extended_defaults(); cfg.plotMode = 'development';
cfg = extended_apply_plot_cfg(cfg);
extended_script10_run(cohort, cfg);   % batch-safe (no questdlg)
run_extended_script11_dominant_periods(cohort)
run_extended_script12_sex_stratified_profiles(cohort)
```

## Script map

| Extended | Replaces Kent | Input | Key outputs |
|----------|---------------|-------|-------------|
| **4** | B ridge + C | Core `Handoff/` | `ExtendedHandoff/AcrossPhotoperiod_Input/WP_*`, `HSubSupported_PeriodMap.mat` |
| **5** | D + D.5 | Script 4 handoff | `Ultradian_RidgePhase_Resync/`, `TransitionEffect_vs_Photoperiod/` |
| **6** | E | Script 4 handoff | LME/FDR workbooks (`AcrossPhotoperiod_LME/`) |
| **7** | F + G | Script 4 handoff | Phase events + profiles (legacy; off in dev by default) |
| **8** | — | Cohort root | `Script8_PublicationFigures_*` composites |
| **9** | — | Cohort root | `Script9_SupplementaryFigures_*` (sex, cluster grids) |
| **10** | — | Cohort root | `Script10_ManuscriptFigures_*` locked Figs 1–7 + `Package_Manuscript/` |
| **11** | — | Cohort root (Script 7 profiles) | `Script11_DominantPeriod_*` mean+median τ + violins |
| **12** | — | Cohort root (Script 7 profiles) | `Script12_SexStratifiedProfiles_*` F\|M activity + post-LD first-peak amplitude |

## Script 11 dominant period outputs

Read-only after Script 7. Uses locked clusters (UR_1_3 C01, UR_3_6 C01+C02) and **`ClusterMembership.RawPeriod_h`** (no `WP_TS`). Per mouse: mean + median of candidate periods; cohort: mean + median across mice.

| Output | Use |
|--------|-----|
| `Tables/DominantPeriod_ClusterSummary.csv` | Cohort median/mean τ + IQR/SD (`IncludeInMainText` flag) |
| `Tables/DominantPeriod_PerMouse.csv` | Per-mouse mean + median audit |
| `Tables/DominantPeriod_Output.xlsx` | Same tables + Settings |
| `Figures/Supp_DominantPeriod_ClusterPeriod_3Cluster.*` | Per-mouse candidate-period violins (magenta median / black mean) |

Script 11 is **standalone** — Script 10 does not copy these into `Package_Manuscript/`.

## Script 12 sex-stratified profiles

Female | Male side-by-side L12–L22 **activity** grids for UR_1_3 C01, UR_3_6 C01, UR_3_6 C02, plus **post–lights-off first-peak amplitude** (F/M overlaid; mean ± SD; Mann–Whitney + BH-FDR stars only on figure). Amplitude = first local max after lights-off minus pre-off baseline (search window ≈ 1× cluster period, capped at 3 h). No coherence / ridge-power panels. Standalone; not in Script 10.

## Script 10 locked main figures

| Fig | Content |
|-----|---------|
| 1 | Circadian characteristics — **placeholder** (collaborator external) |
| 2 | Multiscale RAW + HSub residual **Female\|Male** + CarryForward retention E |
| 3 | Co-expression: A abs power; B descriptive Δ; **C LME forest** (Script 6) |
| 4 | UR_1_3 Z-scored activity L12–L22 (**C01**) |
| 5 | UR_1_3 **24h ZT coherence (A only)** on main |
| 6 | UR_3_6 activity twin of Fig 4 (**C01**) |
| 7 | UR_3_6 coherence twin of Fig 5 (**A only**) |

**Supplemental (Package `Figures/Supp/`, not renumbered mains):**
- `Supp_Transitions/` — Fig05/07 **LD Pre/Post R** bars (Script 5 BH-FDR ΔR>0 stars; cite `Resync_PrimaryStats_BH_FDR` in text)
- `Supp/` — UR_3_6 **C02** activity + coherence twins (`cfg.script10.manuscriptClusters`)

`Package_Manuscript/` includes Figs **3–7** + Supp + **full supplemental tables** (`Tables/Workbooks`, `Tables/CSV`, index); **excludes** Fig 1 and Fig 2. Same table bundle mirrored under `Script10_*/Tables/Supplemental/`. **Script 11** (dominant UR periods + violins) is a **separate** output folder — not bundled into Script 10.

**Do not conflate:** Fig 3C = co-expression LME; Script 5 transition FDR = confirmatory resync (LD Pre/Post R supplemental, not glued onto Figs 5/7). Ridge-power / ΔR gradient figures are not exported.

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
