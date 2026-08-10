# Changelog: BioRhythms_Behaviour

## [Unreleased] — Script 10 manuscript amendments

### Changed
- Figs 5/7 **main composites = panel A only** (24h ZT coherence grid); transition supp = **LD Pre/Post R** bars (`Supp_Fig05/07_LD_PrePost_R`) under `Standalone/Supp_Transitions/` — not glued into Wide/Tall
- Dropped Script 10 export of `Supp_Fig05/07_B_RidgePower` and `_C_DeltaR` gradient figures (Script 5 FDR tables remain in Package for citation)
- UR_3_6 **C02** (ClusterRank 2, ~4.25–5.0 h) activity + coherence twins as Supp extras (`Supp_Activity_UR36_C02`, `Supp_Coherence_UR36_C02`); mains stay C01
- `cfg.script10.manuscriptClusters` lock: resolve by ClusterRank (+ optional period-window warn); builders accept explicit ClusterID
- Package_Manuscript: main Figs 3–7 unchanged set; Supp under `Figures/Supp/`; Script 5 FDR tables still packaged for text citation
- Figs 5/7 panel A: 24 h ZT phase coherence (Script 7 profiles) instead of peri-transition DL/LD curves
- Activity + coherence grids: figure-level `n = … mice (across photoperiods)` plus per-facet `n = …` (cluster membership can be sparse at early L; union ≠ per panel)
- Dark-phase shade + lights-off line aligned with Script 9 / collaborator style

### Fixed
- `cfg.plotMode = 'publication'` now writes `Script10_ManuscriptFigures_Publication` (mode string normalised in runner + `make_output_dirs_`)
- Fig04/Fig06 (and coherence) `n` mice: count unique `SignalID`/`MouseID` across photoperiods — not `File|SignalID` (File differs per L12–L22 and inflated n above study N=24)
- Fig04/Fig06 annotation clarified: bare `n = 24 mice` was a facet union and misleading when early photoperiods plot fewer grey traces (UR_3_6 C01 membership)

## [Unreleased] — Script 7 activity column match (v2.1)

### Fixed
- `find_table_column` in `extended_script7_profiles_run.m`: sex-prefixed SignalIDs (`F_`/`M_`) now strip before matching Excel headers; fuzzy `contains` no longer uses R2025b scalar broadcasting that bound every mouse to `Time (hr)`. Non-activity columns (time/light/ZT) are excluded. Activity z-scored 24h profiles should again show grey individuals + ultradian structure.

## [Unreleased] — Repo layout cleanup

### Changed
- Legacy monolithic scripts moved to `Legacy/` (`harmonic_subtraction_v6.m`, wavelet v4/v5, `PeriodComparison_v4.m`, `Plot_v5.m`)
- Core docs moved into module folders: `Analysis/METHODS.md`, `Functions/io/DATA_STRUCTURE.md`, `Functions/handoff/HANDOFF_SCHEMA.md`, `Functions/hsub/README.md`
- Root retains project-level docs only (`README.md`, `PROJECT_CONTEXT.md`, `CHANGELOG.md`, `TODO.md`, `CURSOR_MATLAB_WORKFLOW.md`, `LabNotebook.md`)
- `setup_paths()` adds `Legacy/` to the MATLAB path

## [Unreleased] — Extended UR E0–E1

### Added
- Top-level `Extended/` tree with Kent A–G copies in `Legacy/Kent_AG/`
- `setup_extended_paths.m`, `Config/extended_defaults.m`, `EXTENDED_CONTEXT.md`, `METHODS_EXTENDED.md`
- Thin E1 wrappers `Analysis/run_extended_*.m` (A–G via `run()`)

### Notes
- Core Scripts 1–3 frozen; CarryForward ±15% SEL_P360; LL projected dark = L22 ZT22–24
- Core Script 3 residual-CWT vs Extended C are parallel claims

## [Unreleased] — Core v1 build

### Added
- Three-script Core v1 pipeline: `run_harmonic_subtraction`, `run_behav_wavelet`, `run_across_photoperiod`
- `Config/core_defaults.m`, `Config/collaborator_palette.m`
- Modular `Functions/{io,hsub,wavelet,handoff,plot}/`
- `HANDOFF_SCHEMA.md`, frozen `PROJECT_CONTEXT.md` / `DATA_STRUCTURE.md` / `METHODS.md`
- `user_paths.example.m` for gitignored results roots

### Changed
- `paths.m` adds `Config/` and recursive `Functions/`; no default results directory
- `README.md` points to Core v1 entry points
- `setup_paths` fprintf notes optional `setup_extended_paths`

### Pending (core-v1.0 tag)
- MATLAB checkpoints Phases 2–6

## [Pre-Core] — Legacy monolithic scripts

- Scripts now in `Legacy/`: `harmonic_subtraction_v6.m`, `behav_wavelet_scalogram_v5.m`, `behav_wavelet_powerspectrum_v4.m`, `PeriodComparison_v4.m`, `Plot_v5.m`
