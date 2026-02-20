# Harmonic_Removal_MiceBehaviour_RegressionAnchor_BlockShuffle_v5

A MATLAB pipeline for detecting and subtracting circadian-scale harmonic structure from mouse behavioural time-series, with QC comparisons and an explicit per-mouse recommendation of which residual signal to take forward.

---

## What this script does

For each selected Data column (each column is treated as an independent mouse):

1. **Reads input data from an Excel file** (preserving original column names) and lets you choose:
   - **Time** column (numeric or datetime)
   - **Light duration condition (h)** column (preserved in outputs)
   - **One or more Data columns**
   - Optional **exclude list** (safe by default, excludes nothing unless you choose)

2. **Handles missing data (conservative policy)**
   - If missingness **≤ 5%**: fills short gaps only using **pchip** interpolation with **MaxGap = 30 min**
   - If missingness **> 5%**: no interpolation, analysis uses available samples

3. **Removes slow baseline drift**
   - Uses a **moving median** baseline with a dynamically chosen window (roughly 72–168 h, scaled to recording length)
   - Performs all modelling on the **detrended** series

4. **Detects a circadian anchor period (22–28 h)**
   - Scans candidate periods in **22–28 h** with **0.01 h step**
   - Compares a **drift-only regression** against **drift + sinusoid**
   - Picks the period with maximum **ΔR²** (DeltaR2)

5. **Gates acceptance of the anchor (`AnchorOK`) using block-shuffle surrogates**
   - Generates a behaviour-friendly null via **block-shuffle** surrogates
   - Computes a surrogate p-value (`pBlock`) from the distribution of surrogate max ΔR²
   - Uses hard viability guards (cycles, ΔR² floor, edge guard), while:
     - `PeakZ` and `AmpSNR` are **QC flags only**, not hard rejection criteria

6. **If `AnchorOK`, runs harmonic subtraction in four modes**
   The script always runs both minimum-period passes:
   - **Min360**: minimum period of interest = 360 min
   - **Min60**: minimum period of interest = 60 min

   For each minimum-period pass, it produces two residuals:
   - **Selective**: removes k=1 and only those higher harmonics k≥2 that behave like true harmonics across windows
   - **FullLadder**: removes the full harmonic ladder k=1..K down to the chosen minimum period

   This yields **4 residual outputs per mouse**:
   - `Residual_FullLadder_Min360`
   - `Residual_FullLadder_Min60`
   - `Residual_Selective_Min360`
   - `Residual_Selective_Min60`

7. **Generates QC comparisons**
   - **Within-min QC**: Selective vs FullLadder (same Min360 or Min60)
     - uses `nRMS` and `DeltaVarExpl_FullMinusSel`
   - **Cross-min QC**: FullLadder Min60 vs Min360
     - uses `nRMS_60vs360` and `DeltaVarExpl_60minus360`
   - Optionally computes a **harmonic susceptibility table** from wavelet global spectrum attenuation (Min60 vs Min360)

8. **Produces a final per-mouse recommendation**
   - Picks one of the 4 residuals per mouse
   - The recommendation is based on:
     - Cross-min QC to choose **Min360 vs Min60**
     - Within-min QC to choose **Selective vs FullLadder**
   - The sheet includes the **metrics and thresholds used** and a human-readable reason string

---

## Input requirements

### Supported file type
- `.xlsx` (Excel)

### Required columns
Your input spreadsheet must contain at least:
- **Time** column (selected via UI)
  - Can be `datetime` or numeric
- **Light duration condition (h)** column (selected via UI, preserved)
- **≥ 1 Data column** (selected via UI)
  - Each Data column is treated as an independent mouse signal

### Time column formats
The script accepts:
- `datetime`
- numeric time vectors where the unit is inferred from the header name:
  - contains `"min"` → minutes
  - contains `"hr"` or `"hour"` → hours
  - contains `"day"` → days
  - otherwise it attempts a heuristic inference

---

## Missing data policy (important)

- **≤ 5% missing**:
  - fills short gaps only using `fillmissing(...,'pchip','MaxGap',maxGapSamples)`
  - `MaxGap` corresponds to **30 minutes**, converted into samples using the inferred sampling interval
- **> 5% missing**:
  - no interpolation
  - all regression and wavelet analyses operate on available samples

Gap handling is logged to the **GapFill_Report** sheet.

---

## Anchor detection and acceptance

### Anchor period search
- Candidate periods: **22–28 h**
- Step: **0.01 h**
- Model comparison:
  - baseline drift model: intercept + linear + quadratic in standardised time
  - augmented model: baseline drift + sinusoid at candidate period
- Score: **ΔR²** (improvement in explained variance)

### Acceptance gate (`AnchorOK`)
An anchor is accepted if all of the following are satisfied:
- block-shuffle p-value: `pBlock < 0.05`
- enough cycles: `duration / P0 ≥ 3.5`
- minimum strength: `ΔR² ≥ 0.005`
- edge guard (enabled): P0 not within 0.10 h of the band edges

QC flags (not hard rejection):
- `PeakZ` (broad peak warning if < 1.0)
- `AmpSNR` (low SNR warning if < 1.0)

If `AnchorOK` is false, subtraction is not performed and outputs are pass-through.

---

## Harmonic subtraction modes

### FullLadder
- Removes all harmonics **k = 1..K**
- `K` is computed from `floor(P0_min / MinPeriod)` and capped to 24

### Selective
- Always removes **k = 1** (fundamental) by default
- Adds higher harmonics k≥2 only if they pass “harmonic-likeness” criteria across windows:
  - period ratio stability (MAD constraint)
  - phase-locking (PLV-style criterion)
  - amplitude coupling (Spearman correlation in long regime or relative amplitude in short regime)

Windowing:
- uses overlapping windows (dynamic length, typically 72–168 h) stepping by 24 h
- the script records window-level fits in the **Detail** workbook

---

## Outputs

### Folder structure created

<OutputFolder>/
	Reports/
		HarmonicRemoval_Reports_Summary.xlsx
		HarmonicRemoval_Reports_Detail.xlsx
	TimeSeries/
		Residual/
			Residual_FullLadder_Min360_<InputFileName>.xlsx
			Residual_FullLadder_Min60_<InputFileName>.xlsx
			Residual_Selective_Min360_<InputFileName>.xlsx
			Residual_Selective_Min60_<InputFileName>.xlsx
		Removed/
			Removed_FullLadder_Min360_<InputFileName>.xlsx
			Removed_FullLadder_Min60_<InputFileName>.xlsx
			Removed_Selective_Min360_<InputFileName>.xlsx
			Removed_Selective_Min60_<InputFileName>.xlsx
	Figures_Anchor/
		Anchor_<Mouse>.jpg
	Figures_Scalograms/
		Residual/
			Scalogram_Res_FullLadder_<MinKey><Mouse>.jpg
			Scalogram_Res_Selective<MinKey><Mouse>.jpg
		Removed/
			Scalogram_Rem_FullLadder<MinKey><Mouse>.jpg
			Scalogram_Rem_Selective<MinKey>_<Mouse>.jpg


### Time-series Excel outputs
Each output workbook contains:
- Time column
- One column per selected Data column, with headers prefixed by the residual/removed mode
- Light duration condition column

Note: the script writes the Light condition column at the end and preserves the original values.

### Figures
- Anchor figures: 600 DPI JPG
- Scalograms (if enabled): 600 DPI JPG
- Harmonics diagnostic figures are intentionally not produced

---

## Reports and how to interpret them

### Summary workbook: `HarmonicRemoval_Reports_Summary.xlsx`

Sheets you will see:

- **README**
  - A run-specific guide to interpretation and thresholds (written by the script)

- **Index**
  - Sheet map

- **Excluded_Columns**
  - Records any Data columns you excluded via the optional exclude UI step

- **GapFill_Report**
  - Missing fraction per mouse
  - Whether interpolation was applied
  - How many samples/gaps were filled and the longest gap

- **Anchor_Report**
  - AnchorOK
  - Period (h), ΔR², PeakZ, amplitude, AmpSNR
  - duration and cycles at the selected period
  - pBlock and notes explaining acceptance or rejection

- **FullLadder_Min360 / FullLadder_Min60**
  - K removed and `VarExplained` for FullLadder per mouse and min-period

- **Selective_Min360 / Selective_Min60**
  - Which k values were removed (as a list)
  - `VarExplained`, number of usable windows, and regime notes

- **QC_SelVsFull_Min360 / QC_SelVsFull_Min60**
  - Within-min comparison of Selective vs FullLadder:
    - `nRMS`: normalised RMS difference between residuals
    - `DeltaVarExpl_FullMinusSel`: extra variance removed by FullLadder vs Selective
    - Risk label (Low/Moderate/High)

- **QC_CrossMin_360vs60**
  - Cross-min comparison (FullLadder Min60 vs Min360):
    - large differences suggest meaningful content in high-order harmonics below 6 h

- **Recommendation**
  - The final residual choice per mouse
  - Provides:
    - the residual label (one of 4)
    - where to find the file
    - the exact metrics and reasoning string

- **Errors**
  - Any caught exceptions per mouse (processing continues for other columns)

### Detail workbook: `HarmonicRemoval_Reports_Detail.xlsx`
- **Sel_WindowLevel_Min360 / Sel_WindowLevel_Min60**
  - Per-window sinusoid fits for fundamental and candidate harmonics
  - Includes best-fit scanned period around each harmonic target

- **Suscept_Periods_360vs60**
  - Wavelet-derived attenuation ratio:
    - `AttenuationRatio = GWS_Min60 / GWS_Min360`
    - `< 0.30` labelled harmonic-susceptible
    - `0.30–0.70` ambiguous
    - `> 0.70` harmonic-robust

---

## Key parameters (defaults in v5)

- Anchor band: 22–28 h
- Anchor step: 0.01 h
- Block shuffle: 8 h blocks, 300 surrogates, alpha = 0.05
- Viability guards: ≥ 3.5 cycles, ΔR² ≥ 0.005, edge margin 0.10 h
- Missingness threshold: 5% (MaxGap 30 min if eligible)
- Minimum periods always run: 360 min and 60 min
- Figures: Anchor = ON, Scalograms = ON
- Image export: JPG, 600 DPI

---

## Practical guidance

### If AnchorOK is often false
- Check recording duration (need enough cycles)
- Check for heavy missingness or long gaps
- Check whether the signal is dominated by non-stationary artefacts that break the sinusoid model

### If Min60 keeps getting recommended
That usually means:
- FullLadder Min60 differs materially from Min360 (large cross-min nRMS or ΔVarExpl)
- Interpretable as high-order harmonic energy below 6 h that affects the residual unless explicitly removed

### If FullLadder is chosen over Selective within a min-period
That typically indicates:
- Selective did not classify enough harmonics as “true” harmonics across windows
- FullLadder removes additional structured content (as shown by large `DeltaVarExpl_FullMinusSel` and/or large nRMS)

---

## Dependencies

This script uses core MATLAB plus Wavelet Toolbox functions if scalograms are enabled:
- `cwtfilterbank`, `cwt`
- `hours`, `minutes`, `datetime`

If scalograms are disabled, the harmonic subtraction and reports still run.

---

## Versioning

- Script name: `Harmonic_Removal_MiceBehaviour_RegressionAnchor_BlockShuffle_v6`
- Outputs are designed to be deterministic per run (workbooks are deleted and recreated at the start of each run)

---

## Isaiah J. Ting. 2026
