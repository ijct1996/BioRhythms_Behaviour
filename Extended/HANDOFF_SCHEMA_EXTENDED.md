# Handoff schema (Extended UR — Kent B → C → D–G)

Extended handoff contract for the **validated-Raw ultradian** pipeline. Core handoff remains in `Functions/handoff/HANDOFF_SCHEMA.md` (`CoreSummary__*.mat`). Do not conflate the two.

## Parallel claim vs Core Script 3

| | Core Script 3 | Extended C (+ D–G) |
|--|---------------|---------------------|
| UR signal | HSub **residual** CWT peaks / band power | **Raw** ridge period candidates |
| Validation | AnchorOK; exclude harmonics of P₀ (±5%, k≤6) | SEL_P360 period match **±15%** → **CarryForward** |
| Downstream | `across_compute_hsub_validated_ur` | `CarryForward_Periods` in `HSubSupported_PeriodMap` |

Report both as separate confirmatory strategies if both appear in a manuscript. See `METHODS_EXTENDED.md`.

---

## Folder: `AcrossPhotoperiod_Input`

Produced by Kent B (`Behav_wavelet_v12` / `run_extended_wavelet_ridge`). Default name under the wavelet output tree:

| Artifact | Pattern | Contents |
|----------|---------|----------|
| Summary MAT | `WP_Summary__{fileStem}.mat` | `pkgS` — tables + meta (small) |
| Time-series MAT | `WP_TS__{fileStem}.mat` | `pkgTS` — long TS tables |
| Sidecar (optional) | `WP_BandColourMap.xlsx` | Band colour map (does not alter pkgS/pkgTS) |

Extended C (`extended_period_gate_run`) reads **only** `WP_Summary__*.mat` from this folder (needs `PeriodCandidates_Long`). WP_TS is consumed by D–G (ridge phase / profiles), not by the period gate.

---

## `WP_Summary__*.mat` — `pkgS`

| Field | Description |
|-------|-------------|
| `pkgS.meta` | Script id, timestamp, Raw/HSub flags, band definitions, candidate QC mins, primary/secondary HSub residual names |
| `pkgS.file` | `FileStem`, `InputFile`, `PerFileFolder` |
| `pkgS.tables.Inputs` | Input inventory |
| `pkgS.tables.Peaks` | Peak summary |
| `pkgS.tables.DetailIndex` | Detail index |
| `pkgS.tables.BandSummary` | Band summary |
| `pkgS.tables.BandConditionSummary` | Band × condition summary |
| `pkgS.tables.PlotMeta` | Scalogram plot metadata |
| **`pkgS.tables.PeriodCandidates_Long`** | **Required by Extended C** — Raw + HSub ridge period candidates |

### `PeriodCandidates_Long` columns (Kent B)

`File`, `SignalID`, `ConditionParsed`, `Source` (`RAW` / `HSUB`), `HSubResidualMode` (e.g. `SEL_P360`, `SEL_P60`, `FL_P360`, `FL_P60`), `Photoperiod_h`, `Phase`, `BandName`, `CandidateID`, `CandidateRank`, `MedianRidgePeriod_h`, `IQR_RidgePeriod_h`, band/ridge power stats, `RidgeCoverageFrac`, `COIValidFrac`, `ValidPointCount`, `TotalPointCount`, `PassQC`, `QCReason`.

---

## `WP_TS__*.mat` — `pkgTS`

| Field | Description |
|-------|-------------|
| `pkgTS.meta` / `pkgTS.file` | Copied from matching `pkgS` |
| `pkgTS.tables.BandPower_Long` | Long band power (+ `ZT_hours`, `ValidFlag`) |
| `pkgTS.tables.RidgePeriod_Long` | Long ridge period |
| `pkgTS.tables.RidgePower_Long` | Long ridge power |
| `pkgTS.tables.RidgeTrajectory_Long` | Ridge trajectory |
| `pkgTS.tables.RidgePhase_Long` | Ridge phase (rad) |

Typical long-row keys: `File`, `SignalID`, `Source`, `BandName`, `Time_days`, `ZT_hours`, `LightStateValue`, `Phase`, `Value` / ridge fields, `ValidFlag`.

---

## Gate output: `HSubSupported_PeriodMap`

Written under:

`{AcrossPhotoperiod_Input}/RawVsSelectiveHSub_PeriodValidation/`

| File | Description |
|------|-------------|
| `HSubSupported_PeriodMap.mat` | `validationMap` struct (`-v7.3`) |
| `HSubSupported_PeriodMap.xlsx` | Same tables as workbook sheets |
| `Figures/` | JPEG validation figures (600 DPI default) |
| `Logs/` | Run log |

### Locked gate defaults (`extended_defaults` / Kent C)

- `PRIMARY_HSUB_MODE` = **SEL_P360**
- `PERIOD_TOLERANCE_FRAC` = **0.15** (±15%; match uses `log2(1 + frac)`)
- Primary phase = **All**
- Ridge coverage / COI ≥ **0.50**; `PassQC` required for Raw and primary HSub unless overridden in cfg

### `validationMap` fields

| Field | Role |
|-------|------|
| `Settings` | Run settings snapshot |
| `LoadSummary` | Per–WP_Summary load status |
| `AllPeriodCandidates` | Concatenated standardised candidates (+ `GlobalCandidateRow`, `SourcePackage`) |
| `Matched_Periods_All` | Every Raw UR candidate with match status / mode columns |
| **`CarryForward_Periods`** | **Downstream D–G inclusion table** (`CarryForward == true`) |
| `RawOnly_NotCarriedForward` | Raw candidates that did not CarryForward |
| `HSubOnly_ResidualFeatures` | Eligible primary HSub without Raw match |
| `FullLadder_Sensitivity` | Slim FL view of matched rows |
| `Retention_ByBand` / `Retention_ByPhotoperiod` / `Retention_ByPhotoperiodBand` | Retention summaries |
| `QC_Flags` | Non-carry / harmonic-sensitive / FL-sensitive rows |

### Key `CarryForward_Periods` / `Matched_Periods_All` columns

Identity: `File`, `SignalID`, `ConditionParsed`, `Photoperiod_h`, `Phase`, `BandName`  
Raw: `RawCandidateID`, `RawPeriod_h`, `RawPassQC`, ridge/COI fields  
Primary HSub: `PrimaryHSubMode`, `PrimaryHSubPeriod_h`, …  
Match: `AbsLog2PeriodDiff`, `PeriodDiffPercent`, `WithinTolerance`, `MatchStatus`, **`CarryForward`**, `HarmonicSensitive12hFlag`, `QCFlag`  
Labels: `FullLadderSensitivityStatus`, `FinalValidationClass`  
Plus per-mode columns for secondary/sensitivity modes (`SEL_P60_*`, `FL_P360_*`, `FL_P60_*`).

---

## Entry points

```matlab
setup_extended_paths
run_extended_period_gate                    % interactive if no args
run_extended_period_gate(handoffDir)        % non-interactive
validationMap = extended_period_gate_run(handoffDir, extended_defaults())
```

Legacy monolith: `Extended/Legacy/Kent_AG/C_Validate_RawVsSelectiveHSub_Periods_v1.m` (immutable reference).

---

## Resync output: `Ultradian_RidgePhase_Resync`

Written under:

`{AcrossPhotoperiod_Input}/Ultradian_RidgePhase_Resync/`

| File | Description |
|------|-------------|
| `Ultradian_RidgePhase_Resync_Output.mat` | `resync` struct (`-v7.3`) |
| `Ultradian_RidgePhase_Resync_Output.xlsx` | Same tables as workbook sheets |
| `Figures/` | Main D figures + `LL_Projected_Aftereffect/` |
| `RidgePower_Takeaway3_Figures/` | Absorbed D.5 ridge-power takeaway JPEG/TIFF panels + summary workbook |
| `Logs/` | Run log |

### Locked D defaults (`extended_defaults` / modular E2)

- Peri-transition window = **6 h**
- Bin width = **0.5 h**
- Pre/post summary window = **2 h**
- Minimum observations per bin / summary = **10**
- Candidate-level permutation tests = **10,000**
- FDR = **Benjamini-Hochberg**, `alpha = 0.05`
- LL projected reference photoperiod = **L22**

### `resync` fields

| Field | Role |
|-------|------|
| `Settings` | Run settings snapshot |
| `LoadSummary` | Per-`WP_TS__*.mat` load status |
| `ValidatedCandidates_Used` | CarryForward candidate lookup keyed by `CandidateID` |
| `TransitionPhase_Long` | Real LD/DL + pseudo-transition long table |
| `BinnedCoherence` / `PrePostCoherence` / `DeltaR_Summary` | Group summaries for real transitions |
| `CandidatePrePost` / `CandidateDeltaR` | Candidate-level pre/post coherence summaries |
| `Resync_PrimaryStats` / `_BH_FDR` | Primary real-transition `DeltaR > 0` family |
| `Resync_RealVsPseudoStats` / `_BH_FDR` | Real-vs-pseudo control family |
| `Resync_RidgePeriodStats` / `_BH_FDR` | Secondary ridge-period pre/post family |
| `Resync_RidgePowerStats` / `_BH_FDR` | Secondary ridge-power pre/post family |
| `LL_NoTransitionSummary` | LL candidates excluded from true transition inference |
| `LL_ProjectedPhase_Long` | Secondary LL projected-event long table |
| `LL_ProjectedBinnedCoherence`, `LL_ProjectedPrePostCoherence`, `LL_ProjectedDeltaR` | LL projected summaries |
| `LL_ProjectedCandidatePrePost`, `LL_ProjectedCandidateDeltaR` | LL projected candidate summaries |
| `LL_ProjectDayDecayStats` / `_BH_FDR` | LL projected day-decay family |
| `L22_vs_LLProjectedStats` / `_BH_FDR` | Final real L22 vs projected LL aftereffect family |
| `outRoot`, `outXLSX`, `outMAT` | Convenience output paths |

### FDR families

Each family is corrected separately within `extended_bh_fdr` using the `FDR_Family` column:

- `Resync_Primary_DeltaR`
- `Resync_RealVsPseudo`
- `Resync_RidgePeriod_PrePost`
- `Resync_RidgePower_PrePost`
- `LL_Projected_Primary_DeltaR`
- `LL_Projected_TransitionVsPseudo`
- `LL_Projected_RidgePeriod_PrePost`
- `LL_Projected_RidgePower_PrePost`
- `LL_Projected_DayDecay`
- `L22_vs_LLProjected_Aftereffect`

### Entry points

```matlab
setup_extended_paths
run_extended_ridge_resync
run_extended_ridge_resync(handoffDirOrParent)
resync = extended_ridge_resync_run(handoffDirOrParent, extended_defaults())
```

Legacy monolith: `Extended/Legacy/Kent_AG/D_Ultradian_RidgePhase_Resync_v4_FDR_LLProjected.m` (immutable reference). D.5 ridge-power plots are absorbed into `Extended/Functions/plot/extended_plot_ridge_power_takeaway.m`.
