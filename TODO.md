# TODO: BioRhythms Core v1

## Phase 0 — Repo hygiene
- [x] `.gitignore` (Data/, Figures/, Output/, `*.mat`, `user_paths.m`)
- [x] Folder skeleton: `Analysis/`, `Functions/{io,hsub,wavelet,handoff,plot}/`, `Config/`, `Legacy/`
- [x] `paths.m` updated (code paths only)
- [x] `user_paths.example.m`
- [x] `CURSOR_MATLAB_WORKFLOW.md` → Dev path

## Phase 1 — Documentation freeze
- [x] `PROJECT_CONTEXT.md`
- [x] `DATA_STRUCTURE.md`
- [x] `METHODS.md`
- [x] `HANDOFF_SCHEMA.md`
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
- [ ] Move v4–v6 scripts to `Legacy/`
- [x] `Legacy/README.md`
- [x] Root `README.md` → three `run_*.m` entry points
- [ ] `CHANGELOG.md` → core-v1.0
