# BioRhythms_Behaviour — Core v1 + Extended UR

Mouse locomotor rhythms: harmonic subtraction (validation), RAW wavelet, across-photoperiod co-expression. Extended UR (Kent A–G) is a parallel validated-Raw ridge pipeline under `Extended/`.

## Quick start — Core v1 (MATLAB R2025b)

```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour')
setup_paths();
copyfile('user_paths.example.m', 'user_paths.m');  % once — edit outputRoot

run_harmonic_subtraction   % Script 1
run_behav_wavelet          % Script 2
run_across_photoperiod     % Script 3
```

## Quick start — Extended UR (E1)

```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour\Extended')
setup_extended_paths   % also pulls in Core setup_paths

run_extended_hsub                 % A
run_extended_wavelet_ridge        % B
run_extended_period_gate          % C
run_extended_ridge_resync         % D
run_extended_across_photoperiod   % E
run_extended_phase_events         % F
run_extended_publication_profiles % G
```

See `Extended/README.md`. Core Scripts 1–3 are frozen; do not conflate Core Script 3 residual-CWT UR with Extended CarryForward (±15% SEL_P360).

## Folder layout

```
BioRhythms_Behaviour/
├── README.md, PROJECT_CONTEXT.md, CHANGELOG.md, TODO.md, …   ← project-level docs
├── setup_paths.m, paths.m, user_paths.example.m              ← entry points
├── Analysis/          run_*.m scripts + METHODS.md
├── Config/            core_defaults, palettes
├── Functions/         hsub/, wavelet/, handoff/, plot/, io/ (+ module README / schema)
├── Legacy/            pre–Core v1 monolithic scripts
└── Extended/          UR ridge pipeline (Scripts 4–8) + Extended docs
```

## Docs

| File | Purpose |
|------|---------|
| `PROJECT_CONTEXT.md` | Science + scope (Core + Extended) |
| `Analysis/METHODS.md` | Core pipeline logic |
| `Extended/METHODS_EXTENDED.md` | Extended CarryForward / projected dark |
| `Functions/io/DATA_STRUCTURE.md` | Pipeline_Input layout |
| `Functions/handoff/HANDOFF_SCHEMA.md` | Script 2 → 3 contract |
| `Functions/hsub/README.md` | Harmonic subtraction deep dive |
| `CURSOR_MATLAB_WORKFLOW.md` | Cursor + MATLAB loop |

## Data (ready)

`C:\Users\User\Dev\Cursor\Research\Chronobiology\Data\BioRhythms Behaviour Data\Pipeline_Input\`

## Legacy

Pre–Core v1 scripts live in `Legacy/` until Phase 7 regression passes — see `Legacy/README.md`. Kent A–G copies: `Extended/Legacy/Kent_AG/`.
