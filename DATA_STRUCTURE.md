# Data Structure: BioRhythms Core v1

> Raw data **not** in git. Pipeline_Input files are pre-split and ready — do not re-split unless sources change.

## Storage locations

| Dataset | Path | Git |
|---------|------|-----|
| Pipeline_Input | `C:\Users\User\Dev\Cursor\Research\Chronobiology\Data\BioRhythms Behaviour Data\Pipeline_Input\` | No |
| Manifest | `...\PIPELINE_INPUT_MANIFEST.csv` | No |
| Kent originals | `C:\Users\User\OneDrive\Desktop\Files\Kent\04. Data\03. Behaviour_Mice\LB_Thesis_Data_LP\` | No |
| Results | User-selected (see `user_paths.m`) | No |

## Cohorts (Core v1)

| Folder | Files | Rows | Mice | Duration |
|--------|-------|------|------|----------|
| `C57/` | L12–L24 (7) | 721 | 25 | ~5 d |
| `NR2B_LP/` | L12–L24 (7) | 721 | 15 | ~5 d |
| `NR2B_LD_DD/` | L12_NR2B_LD, L0_NR2B_DD | 1441 | 36 | ~10 d |

## Filename convention

```
L{photoperiod}_{genotype}[_LP|_LD|_DD]_02JUL26.xlsx
```

Examples: `L12_C57_02JUL26.xlsx`, `L24_NR2B_LP_02JUL26.xlsx`, `L0_NR2B_DD_02JUL26.xlsx`

## Column layout (every file)

| Order | Column | Type | Notes |
|-------|--------|------|-------|
| 1 | `Time (hr)` | numeric | Hours from start; 10 min steps |
| 2…n−1 | Mouse IDs | numeric | One column per animal |
| last | `Light duration (h)` | numeric | Constant per file (photoperiod label) |

## Group definitions

| Cohort | Groups | Assignment |
|--------|--------|------------|
| C57_LP | C57 (single pool) | Auto |
| NR2B_LP | Male / Female if parseable | Column name suffix |
| NR2B_LD_DD | Male / Female if parseable | Column name suffix |

Fallback: interactive group UI (`group_assignment_dialog`).

## Re-split (only if sources change)

Script: `split_pipeline_input.m` (in Data folder). Manifest: `PIPELINE_INPUT_MANIFEST.csv`.

## Handoff outputs (Script 2 → 3)

See `HANDOFF_SCHEMA.md`. No full CWT arrays in MAT files.
