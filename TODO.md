# TODO: BioRhythms Core v1 + Extended UR

## Phase 0 — Repo hygiene
- [x] `.gitignore` (Data/, Figures/, Output/, `*.mat`, `user_paths.m`)
- [x] Folder skeleton: `Analysis/`, `Functions/{io,hsub,wavelet,handoff,plot}/`, `Config/`, `Legacy/`
- [x] `paths.m` updated (code paths only)
- [x] `user_paths.example.m`
- [x] `CURSOR_MATLAB_WORKFLOW.md` → Dev path

## Phase 1 — Documentation freeze
- [x] `PROJECT_CONTEXT.md`
- [x] `Functions/io/DATA_STRUCTURE.md`
- [x] `Analysis/METHODS.md`
- [x] `Functions/handoff/HANDOFF_SCHEMA.md`
- [x] `TODO.md` (this file)

## Phase 2 — Shared foundation
- [x] `Config/core_defaults.m`
- [x] `Config/collaborator_palette.m`
- [x] `Functions/io/*`
- [ ] **Checkpoint**: metadata parse on `L12_C57` and `L0_NR2B_DD` in MATLAB

## Phase 3 — Script 1
- [x] `Analysis/run_harmonic_subtraction.m` + `Functions/hsub/*`
- [ ] **Checkpoint**: `L12_C57_02JUL26.xlsx` → sensible AnchorOK + SEL_P360 recommendation

## Phase 4 — Script 2
- [x] `Analysis/run_behav_wavelet.m` + `Functions/wavelet/*` + handoff
- [ ] **Checkpoint**: L12 C57 → scalogram + `CoreSummary__L12_C57_02JUL26.mat`

## Phase 5 — Script 3
- [x] `Analysis/run_across_photoperiod.m`
- [ ] **Checkpoint**: C57_LP handoff (7 photoperiods) → across-photoperiod outputs

## Phase 6 — Full cohort regression
- [ ] Script 1+2: all C57_LP (7), NR2B_LP (7), NR2B_LD_DD (2)
- [ ] Script 3: three cohorts separately
- [ ] Pause OneDrive during batch

## Phase 7 — Legacy cleanup
- [x] Move v4–v6 scripts to `Legacy/`
- [x] `Legacy/README.md`
- [x] Root `README.md` → three `run_*.m` entry points
- [ ] `CHANGELOG.md` → core-v1.0

---

## Extended UR — E0 / E2 (Kent A–G)

Core Scripts 1–3 remain **frozen**. Do not modify Core `run_harmonic_subtraction` / `run_behav_wavelet` / `run_across_photoperiod` logic for Extended.

### E0 — Scaffold
- [x] Top-level `Extended/` tree (`Analysis/`, `Functions/{handoff,ridge,stats,plot}/`, `Config/`, `Legacy/Kent_AG/`)
- [x] Copy Kent Final Versions A–G into `Extended/Legacy/Kent_AG/` (`D.5` → `D5_Plotting_Ridge_Periods.m`)
- [x] `Extended/Legacy/Kent_AG/README.md` (A–G → wrapper map)
- [x] `Extended/EXTENDED_CONTEXT.md`, `Extended/METHODS_EXTENDED.md`
- [x] `Extended/Config/extended_defaults.m`, `Extended/README.md`
- [x] `Extended/setup_extended_paths.m`

### E1 — Thin wrappers
- [x] `run_extended_hsub` → A
- [x] `run_extended_wavelet_ridge` → B
- [x] `run_extended_period_gate` → C (±15% SEL_P360 CarryForward)
- [x] `run_extended_ridge_resync` → D
- [x] `run_extended_across_photoperiod` → E
- [x] `run_extended_phase_events` → F
- [x] `run_extended_publication_profiles` → G
- [ ] **Checkpoint**: MATLAB smoke — `setup_extended_paths` then one wrapper opens Kent dialogs

### E2 — Modularise
- [x] C CarryForward gate extracted into `Extended/Functions/handoff/`
- [x] D ridge-phase resync extracted into `Extended/Functions/ridge/extended_ridge_resync_run.m`
- [x] BH/FDR helper extracted into `Extended/Functions/stats/extended_bh_fdr.m`
- [x] LL projected event definitions extracted into `Extended/Functions/ridge/extended_ll_projected_event_definitions.m`
- [x] D.5 ridge-power takeaway absorbed into `Extended/Functions/plot/extended_plot_ridge_power_takeaway.m`
- [x] `run_extended_period_gate` and `run_extended_ridge_resync` prefer modular entries with Kent fallback
- [ ] **Checkpoint**: MATLAB smoke — `run_extended_period_gate` then `run_extended_ridge_resync` on a real Extended handoff

### Scripts 4–7 (streamlined Extended pipeline)
- [x] `run_extended_script4_ridge_validation` — ridge handoff from Core Handoff + CarryForward gate
- [x] `run_extended_script5_transition_resync` — transition resync, FDR, photoperiod gradient, stable-days sensitivity
- [x] `run_extended_script6_across_lme` — Kent E legacy (LME/FDR tables)
- [x] `run_extended_script7_phase_profiles` — Kent F+G merged (skipped in dev unless `generateLegacyFigures=true`)
- [x] Development plot mode default (96 DPI PNG); publication deferred
- [x] Primary bands locked: `UR_1_3`, `UR_3_6`
- [ ] **Checkpoint**: MATLAB smoke — Scripts 4→5 on real Core Handoff cohort

### E3 — Remaining modular ports
- [ ] Modularise Kent E (replace legacy 600 DPI runner)
- [ ] Modularise Kent F+G with development/publication plot mode
