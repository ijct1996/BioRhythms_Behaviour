# Changelog: BioRhythms_Behaviour

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
