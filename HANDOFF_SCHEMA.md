# Handoff schema (Core v1)

`handoff_version = '1.0'` — Script 2 → Script 3 contract.

## Files

| File | Location | Description |
|------|----------|-------------|
| `CoreSummary__{fileStem}.mat` | `{cohortRoot}/Handoff/` | Per-photoperiod summary |
| `CoreHandoff_Index.xlsx` | Same | Auto-appended index when output root reused |

## MAT struct fields (`CoreSummary__*.mat`)

| Field | Type | Description |
|-------|------|-------------|
| `handoff_version` | char | `'1.0'` |
| `fileStem` | char | e.g. `L12_C57_02JUL26` |
| `created` | datetime | Write timestamp |
| `meta` | struct | From `parse_file_metadata` + column indices |
| `groups` | struct array | `.name`, `.colIdx` |
| `hsub` | struct | `.available`, `.path`, `.anchor`, `.recommendation` tables |
| `periodTable` | table | Peaks: SignalID, Group, PeakRank, PeakPeriod_hr, PeakValue_log10, LightDuration_h |
| `bandPower` | table | Per mouse: SignalID, Group, LightDuration_h, CR_20_28, UR_1_3, … |
| `figurePaths` | cellstr | Scalogram file paths |
| `paths` | struct | `.raw`, `.hsub`, `.wavelet` |
| `analysisNote` | char | Documents RAW-primary rule |

**Excluded**: full CWT arrays (`wt`, `powerSpec` volumes).

## Index columns (`CoreHandoff_Index.xlsx`)

| Column | Description |
|--------|-------------|
| FileStem | `L12_C57_02JUL26` |
| LightDuration_h | Photoperiod numeric |
| Cohort | `C57_LP`, `NR2B_LP`, `NR2B_LD_DD` |
| SummaryPath | Full path to MAT |
| Created | datetime |

## Co-expression derived fields (Script 3)

Computed in `across_compute_coexpression` from `bandPower`:

- `CR_to_UR_1_3`, `CR_to_UR_3_6`, `CR_to_UR_6_12`, `CR_to_UR_12_18`
- `CR_to_UR_total`
- `UR_*_fraction_of_UR_total`

## Script 3 stitched scalograms

`across_plot_stitched_scalograms` reads `paths.raw` and `groups` from each summary (not stored CWT). Group-mean RAW signals are concatenated in ascending `LightDuration_h` order; photoperiod boundaries are drawn at segment joins.

## Script 3 HSub-validated ultradian (confirmatory)

Computed in `across_compute_hsub_validated_ur` — **does not modify Script 1 or 2**.

| Source | Metric |
|--------|--------|
| RAW xlsx (`paths.raw`) | CR₂₀₋₂₈ band power |
| SEL_P360 residual xlsx (`paths.hsub/{fileStem}/TimeSeries/Residual/`) | UR band power, UR period peaks |
| `hsub.anchor` in handoff | AnchorOK gate, P₀ for harmonic rejection |

Harmonic rule: exclude peak period P if |P − P₀/k| / (P₀/k) ≤ 0.05 for any k = 1..6 within UR range (1–18 h).
