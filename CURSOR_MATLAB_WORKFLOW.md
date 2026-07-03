# MATLAB + Cursor Workflow: BioRhythms_Behaviour

This guide replaces the old ChatGPT copy-paste loop with a file-based, reproducible workflow.

## Old workflow (ChatGPT)

```
MATLAB open → copy code → paste in ChatGPT → copy answer → paste in MATLAB → run
Error? → copy error from Command Window → paste in ChatGPT → repeat
```

**Problems:** no version history, no project memory, context lost each session, easy to paste wrong version.

## New workflow (Cursor + MATLAB + GitKraken)

```
┌─────────────┐     ┌──────────────┐     ┌─────────────┐
│   Cursor    │     │    MATLAB    │     │  GitKraken  │
│ edit .m +   │────▶│ run & test   │────▶│ commit/push │
│ docs        │     │ see figures  │     │ to GitHub   │
└─────────────┘     └──────────────┘     └─────────────┘
       ▲
       │ reads PROJECT_CONTEXT.md, METHODS.md (no re-prompting)
```

**Cursor edits files directly** — you do not copy-paste code into chat for every change.

---

## Session setup (once per day, ~2 minutes)

### 1. Open the project in Cursor
- Open folder: `Research\Chronobiology\BioRhythms_Behaviour`
- Or open `Research\Research.code-workspace`

### 2. In MATLAB
```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour')
setup_paths();   % adds project + Shared/MATLAB to path
theme = plot_config('development');
```

### 3. Read context (you or Cursor)
- `PROJECT_CONTEXT.md` — what this project is about
- `TODO.md` — today's priority
- `METHODS.md` — pipeline logic

---

## How to ask Cursor for help (instead of copy-paste)

### Refactor or improve code
> Read `PROJECT_CONTEXT.md` and `harmonic_subtraction_v6.m`. Refactor the harmonic gating section into a separate function without changing behaviour.

Cursor opens and edits the file **in place**. You run it in MATLAB.

### Fix an error
**Old way:** copy entire error from Command Window → paste in ChatGPT.

**New way:** paste the error in Cursor chat **and reference the file**:
> `@harmonic_subtraction_v6.m` — I get this error when running line 240:  
> `Error using ...`  
> Fix without changing the block-shuffle logic.

Or select the relevant lines in the editor and use **Ctrl+L** (chat with selection).

### Add a feature
> Read `METHODS.md`. Add export of the residual signal to `Output/` after harmonic subtraction. Use `PLOT_MODE = 'development'` for any new figures.

### Review statistics
> Review `PeriodComparison_v4.m` against `Shared/Research_Knowledge.md`. Flag any multiple-comparison issues.

---

## Typical work loop

| Step | Where | Action |
|------|-------|--------|
| 1 | Cursor | Update `TODO.md` with today's task |
| 2 | Cursor | Edit one script (one logical change) |
| 3 | MATLAB | Run script; check figures (development quality, fast) |
| 4 | Cursor | If error, paste error + `@filename` — fix in file |
| 5 | Cursor | Update `CHANGELOG.md` or `LabNotebook.md` if important |
| 6 | GitKraken | Commit with intent message → push |

**One logical change per commit.** Examples:
- `Improve harmonic gating QC flags in harmonic_subtraction_v6`
- `Fix y-axis scaling in behav_wavelet_scalogram_v5`

---

## Plotting during development vs publication

```matlab
% While coding (default):
PLOT_MODE = 'development';
theme = plot_config(PLOT_MODE);
apply_plot_theme(gca, theme);
export_figure(gcf, 'preview_scalogram', theme);   % 96 DPI, fast

% When manuscript-ready:
theme = plot_config('publication');
export_figure(gcf, 'Figures/export/Fig2_scalogram.png', theme);  % Tol bright, TNR, 600 DPI
```

---

## What lives where

| File | Purpose |
|------|---------|
| `PROJECT_CONTEXT.md` | Scientific context — Cursor reads this first |
| `METHODS.md` | Pipeline logic (not code) |
| `DATA_STRUCTURE.md` | Excel columns, groups, units |
| `TODO.md` | Current priorities |
| `HarmonicSubtraction_README.md` | Deep dive on harmonic subtraction |
| `Shared/Research_Knowledge.md` | Lab-wide stats/plotting rules |

You **do not** need to re-explain the project each session if these files are kept up to date.

---

## MATLAB + Cursor side-by-side layout

1. **Monitor 1:** Cursor (code + chat)
2. **Monitor 2:** MATLAB (Command Window + figures)

Or split screen: Cursor left, MATLAB right.

**Tip:** Save in Cursor before running in MATLAB (`Ctrl+S`). MATLAB reads from disk.

---

## When to use GitKraken

- After a working change: **commit + push**
- End of session: always push (backup to GitHub)
- Never use GitHub Desktop day-to-day — GitKraken only

---

## Core v1 entry points

```matlab
run_harmonic_subtraction   % Script 1 — 01_HarmonicSubtraction/
run_behav_wavelet          % Script 2 — 02_Wavelet/ + Handoff/
run_across_photoperiod     % Script 3 — reads Handoff/ only
```

Copy `user_paths.example.m` → `user_paths.m` and set `outputRoot` (e.g. `C:\ResearchData\Results\Core`).

Data: `...\Chronobiology\Data\BioRhythms Behaviour Data\Pipeline_Input\`
