# Changelog: BioRhythms_Behaviour

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

### Pending (core-v1.0 tag)
- MATLAB checkpoints Phases 2–6
- Move legacy v4–v6 to `Legacy/`

## [Pre-Core] — Legacy monolithic scripts

- `harmonic_subtraction_v6.m`, `behav_wavelet_scalogram_v5.m`, `behav_wavelet_powerspectrum_v4.m`, `PeriodComparison_v4.m`, `Plot_v5.m`
