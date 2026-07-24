# Behavioural rhythm analysis — handover package

Self-contained MATLAB codes and methods notes for:

1. Circadian harmonic subtraction  
2. Continuous wavelet transform of locomotor activity  
3. Across-photoperiod co-expression and period summaries  

Prepared so you can run the analysis and pull material for a thesis Methods section. GB English.

**MATLAB target:** R2025b (Wavelet Toolbox required for CWT / `cwtfilterbank`).

---

## What’s in this folder

| File | Role |
|------|------|
| `LB_01_HarmonicSubtraction.m` | Harmonic subtraction only |
| `LB_02_Wavelet_and_AcrossPhotoperiod.m` | Wavelet (Part A) + across-photoperiod summaries (Part B) |
| `README.md` | This file — how to run, inputs/outputs, mechanics |
| `LB_Methods_Information_for_Thesis_v3.docx` | Methods **information** (Desktop prose + full handover deliverables) — preferred thesis source material |
| `Template_Input/` | Example activity file + column-layout notes (`L12_TEMPLATE_example.xlsx`) |

Both `.m` files keep their helper functions at the **bottom of the same file**. You do not need separate function files on the path.

---

## Template input

See `Template_Input/README.md`. Quick start: open `Template_Input/L12_TEMPLATE_example.xlsx` with `LB_01` (synthetic L12 data, 5 days, four mice). Real analyses use the same column order: `Time (hr)`, mouse columns, `Light duration (h)`.

---

## Recommended run order

1. Prepare activity Excel files (one photoperiod / light-duration condition per file).
2. Run `LB_01_HarmonicSubtraction.m` on the cohort.
3. Run `LB_02_…` **Part A** (wavelet) on the same raw files, pointing at the harmonic-subtraction output root.
4. After all photoperiods have an `AnalysisSummary__*.mat` in `Summaries/`, run **Part B** (across photoperiod).

Do not run Part B before Part A is finished for every photoperiod you intend to include.

---

## Input file format

Each workbook should contain at least:

- A **time** column (hours preferred; e.g. `Time (hr)`)
- A **light duration** column (hours; e.g. `Light duration (h)`) — constant within a file for a given photoperiod
- One or more **mouse activity** columns (numeric counts / activity)

Sampling is normally **10 min**. The codes infer the step from the time column when possible.

A worked example ships in `Template_Input/L12_TEMPLATE_example.xlsx`.

---

## Harmonic subtraction (`LB_01_…`)

### What it does

For each mouse column:

1. Light cleanup / optional short-gap fill, then slow baseline removal (moving median).
2. **Circadian anchor scan** over 22–28 h: find the period that maximises ΔR² for a sinusoid-plus-drift model vs drift alone.
3. **Surrogate test**: block-shuffle the detrended series (default 8 h blocks, 300 surrogates). Accept the anchor if p &lt; 0.05, enough cycles (≥ 3.5), ΔR² above a small floor, and the peak is not stuck on the band edge.
4. If the anchor fails → **pass-through** (raw kept; recommendation `PASS`).
5. If the anchor passes → fit sliding windows, screen harmonics of P₀, then subtract:
   - **Selective**: fundamental (+ harmonics that pass the likelihood rules)
   - **Full ladder**: fundamental through K ≈ floor(P₀ / Pmin)
6. Two period floors: **Min360** (6 h) and **Min60** (1 h). Default recommendation when AnchorOK: **SEL_P360** (selective, 360 min floor).

### Outputs (per input file stem)

Under `HarmonicSubtraction/{fileStem}/`:

```
Reports/HarmonicRemoval_Summary.xlsx
  sheets: Anchor_Report, Recommendation
TimeSeries/Residual/
  Residual_Selective_Min360_{stem}.xlsx   ← used by Part B confirmatory layer
  Residual_Selective_Min60_{stem}.xlsx
  Residual_FullLadder_Min360_{stem}.xlsx
  Residual_FullLadder_Min60_{stem}.xlsx
  Residual_Recommended_{stem}.xlsx
TimeSeries/Removed/
  Removed_Selective_*.xlsx / Removed_FullLadder_*.xlsx
Figures_Scalograms/Group_*/                 ← SEL_P360 Removed|Residual (two-panel)
  HSub_Removed_Residual_SEL_P360_{group}_L*.*
  Individuals/…                             ← optional (dialog)
```

Residual columns are named `Residual_Selective_<MouseID>` (and likewise for other arms).

On the first file you also define **condition groups** (name + columns) for the scalogram averages; later files reuse those mouse names. Wavelet Toolbox is required for `Figures_Scalograms/`.

### Key parameters

| Parameter | Value |
|-----------|--------|
| Anchor band | 22–28 h (0.01 h step) |
| Surrogates | 300 block-shuffles, 8 h blocks |
| α (anchor) | 0.05 |
| Min cycles at P₀ | 3.5 |
| Default residual | SEL_P360 |

Set `rng(1)` at the top of the code (already present) before any publication re-run so surrogate draws are reproducible.

---

## Wavelet (`LB_02_…` Part A)

### What it does

- Continuous wavelet transform (**Morlet / `amor`**) on **RAW** activity.
- Period limits **60–1590 min**.
- Group-average and individual **jet** scalograms (RAW).
- **HSub Removed | Residual** two-panel jet scalograms from selective Min360 outputs of `LB_01` (group averages and individuals; **AnchorOK mice only** for averages).
- Top **5** peaks from the mean log₁₀ power spectrum per RAW series.
- Band-mean power for co-expression: CR₂₀₋₂₈, UR₁₋₃, UR₃₋₆, UR₆₋₁₂, UR₁₂₋₁₈.
- Reads `HarmonicRemoval_Summary.xlsx` when present and stores Anchor/Recommendation tables in the summary MAT.

### Outputs

```
Wavelet/{fileStem}/
  Scalograms/Group_*/...                          ← RAW
  Scalograms/HSub_SEL_P360/Group_*/...            ← Removed | Residual (two-panel)
  PeriodPeaks/PeriodPeaks_{stem}.xlsx             ← RAW peaks
Summaries/AnalysisSummary__{stem}.mat             ← required input for Part B
Summaries/Summaries_Index.xlsx                    ← cohort index of summaries
```

`Summaries/` is created beside `Wavelet` (same parent folder).

### Condition groups

On the first file you are asked:

1. How many condition groups  
2. For each group: the group name, then a list dialog to tick which mouse columns belong to it  

Later files reuse the same mouse **names** where possible.

---

## Across photoperiod (`LB_02_…` Part B)

### What it does

Loads every `AnalysisSummary__*.mat` in the chosen `Summaries/` folder (sorted by light duration).

**RAW layer (exploratory)**  
Uses Part A band power and period peaks as recorded on RAW.

**Validated layer (confirmatory)**  
For **AnchorOK** mice only:

- Re-runs CWT on the selective Min360 residual (`Residual_Selective_Min360_*.xlsx`)
- **UR** band power and UR period peaks (1–18 h) from that residual
- **CR₂₀₋₂₈** band power still taken from **RAW** CWT
- Drops period candidates within ±5% of P₀/k for k = 1…6 (harmonic filter)
- Writes a validation log of which candidates were flagged as harmonics

### Outputs

```
AcrossPhotoperiod_{cohort}/
  Tables/
    PeriodPeaks_all_photoperiods_RAW.xlsx
    Coexpression_ratios_RAW.xlsx
    BandPower_hsub_validated.xlsx
    PeriodPeaks_hsub_validated.xlsx
    PeriodPeaks_hsub_validation_log.xlsx
    Coexpression_ratios_hsub_validated.xlsx
  Figures/
    Coexpression_CR_UR_RAW.*
    Period_comparison_RAW.*
    Coexpression_CR_UR_hsub_validated.*
    Period_comparison_hsub_validated.*
    Period_Distributions/
      Period_distribution_RAW_{group}.*
      Period_distribution_hsub_validated_{group}.*
    Stitched_Scalograms/
      Stitched_Scalogram_Average_{group}_{cohort}.*
    Stitched_Scalograms_HSub/
      Stitched_HSub_Removed_Residual_SEL_P360_{group}_{cohort}.*
```

Development export: PNG ~96 DPI. Publication: JPEG ~600 DPI.

---

## Mechanics in one picture

```
Activity Excel files (one photoperiod each)
        │
        ▼
 LB_01  Harmonic subtraction
        → Reports + TimeSeries (Residual/Removed) + Figures_Scalograms/
        │
        ▼
 LB_02 Part A  Wavelet on RAW + HSub Removed|Residual scalograms
               → Wavelet/… + Summaries/AnalysisSummary__*.mat
        │
        ▼
 LB_02 Part B  Across photoperiod
        ├── Tables (RAW + residual-validated)
        ├── Co-expression / period comparison
        ├── Period_Distributions/
        └── Stitched_Scalograms/ (+ _HSub)
```

---

## Dependencies / notes

- **Wavelet Toolbox** (`cwt`, `cwtfilterbank`) — required for `LB_01` scalograms and all of `LB_02`.
- Pause OneDrive sync on the results folder during long batches if files are being locked mid-write.
- If a mouse fails the circadian anchor, the confirmatory UR layer simply **omits** that mouse; RAW exploratory summaries can still include them from Part A.
- Folder names use neutral handover labels (`HarmonicSubtraction/`, `Summaries/`, `AcrossPhotoperiod_*`); deliverables match the full analysis set described in this README and in `LB_Methods_Information_for_Thesis_v3.docx`.
- For thesis wording and equations, use `LB_Methods_Information_for_Thesis_v3.docx` — adapt into your own Methods voice; the document is source material, not a paste-ready Methods section.
- To regenerate the v3 Word file after edits: run `_build_methods_v3_ooxml.ps1` (no Python required).
