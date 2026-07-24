# Template input dataset

This folder shows the Excel/CSV layout expected by `LB_01_HarmonicSubtraction.m` and `LB_02_Wavelet_and_AcrossPhotoperiod.m`.

## Required columns (order)

| Order | Column | Type | Notes |
|-------|--------|------|-------|
| 1 | `Time (hr)` | numeric | Hours from recording start; normally **10 min** steps (0, 0.1667, 0.3333, …) |
| 2 … n−1 | Mouse IDs | numeric | One column per animal (e.g. `F01`, `M02`) |
| last | `Light duration (h)` | numeric | **Constant** within the file — photoperiod label (e.g. 12 for L12) |

One photoperiod (light duration) per file. Do not mix L12 and L16 in the same workbook.

## Files here

| File | What it is |
|------|------------|
| `L12_TEMPLATE_example.xlsx` | Synthetic L12 example (5 days, four mice) for a first dry run |
| `column_layout_example.csv` | Same column layout as a plain CSV |

## Quick start

1. Run `LB_01_HarmonicSubtraction.m` on `L12_TEMPLATE_example.xlsx`.
2. Run `LB_02_…` Part A on the same file, pointing at the harmonic-subtraction output root.
3. For a full across-photoperiod run you need several photoperiod files (each with its own light duration) and Part B after all Part A summaries exist.
