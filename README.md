# BioRhythms_Behaviour — Core v1

Mouse locomotor rhythms: harmonic subtraction (validation), RAW wavelet, across-photoperiod co-expression.

## Quick start (MATLAB R2025b)

```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour')
setup_paths();
copyfile('user_paths.example.m', 'user_paths.m');  % once — edit outputRoot

run_harmonic_subtraction   % Script 1
run_behav_wavelet          % Script 2
run_across_photoperiod     % Script 3
```

## Docs

| File | Purpose |
|------|---------|
| `PROJECT_CONTEXT.md` | Science + scope |
| `METHODS.md` | Pipeline logic |
| `DATA_STRUCTURE.md` | Pipeline_Input layout |
| `HANDOFF_SCHEMA.md` | Script 2 → 3 contract |
| `CURSOR_MATLAB_WORKFLOW.md` | Cursor + MATLAB loop |

## Data (ready)

`C:\Users\User\Dev\Cursor\Research\Chronobiology\Data\BioRhythms Behaviour Data\Pipeline_Input\`

## Legacy

Pre–Core v1 scripts (`harmonic_subtraction_v6.m`, etc.) remain at repo root until Phase 7 regression passes — see `Legacy/README.md`.
