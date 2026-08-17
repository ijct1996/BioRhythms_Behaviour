# Changelog: BioRhythms_Behaviour

## [Unreleased] — Script 14 publication stitched scalograms (v1.0)

### Added
- **Script 14:** Recomputes sex-split stitched RAW + HSub SEL_P360 residual scalograms from Handoff/
- Unified 1280×640 layout, Arial, photoperiod labels above heatmap, sex-pooled clim per signal class
- Outputs under `Script14_Scalograms_{Publication|Development}/`

## [Unreleased] — Script 10/11 publication figure pass

### Changed
- **Script 10** (`plotMode = 'publication'`): Arial font; JPEG @ 600 DPI only (no dual PDF export)
- Fig03: titles/line-end tags/n annotations removed; A y-label `Mean band power (au)`; B `\Delta (UR - CR) log10 power`; C y-label `Variables`, `\Delta` tick labels, BH footnote removed (stars kept)
- Fig02 E: title removed; y `Validated Ultradian Periods (%)` [0–100]; x-label `Ultradian Bands` (tick labels unchanged); scalogram standalones drop MATLAB `title()` in publication mode
- Activity/coherence grids (Figs 4–7): no `sgtitle`/union n; per-panel `n` top-left; x `Zeitgeber Time (ZT; h)` 0–24 (ticks 0:4:24) on all panels; activity y `[-3, 3]` with clip warning; coherence y `[0, 1]`
- **Script 11 v2.5**: publication titles stripped; y-label `Ultradian Period (h)`, ylim `[0, 6]`; Arial in publication mode

## [Unreleased] — Script 11 violins only if n≥6 (v2.2)

### Changed
- Script 11: KDE violins only when n≥6; n≤5 shows points + mean/median ticks (no density for tiny n)
- Removed soft-clip attenuation that made jagged/exaggerated tails; smooth KDE with bw capped to data span
- y-limits still from observed periods only (KDE support does not inflate ylim)

## [Unreleased] — Script 11 pointed violin tips restored

### Changed
- Violin KDE support again extends ~3× bandwidth past the data so tips taper (not flat-trimmed)
- Mild bandwidth cap retained; n≥6 rule unchanged

## [Unreleased] — Script 11 violin bandwidth + N clarification

### Changed
- Violins: bandwidth capped to data span; axis limits from observed periods only (no empty long tails)
- `ByCluster`: added `N_Candidates`, `N_Mice_Activity` (ActivityComponent unique mice; may be < period `N`)
- `ByMouse`: added `N_Candidates` (period values; identical periods can overlap visually)
- Documented: Script 11 `N` ≠ activity-plot panel n (photoperiod-scoped activity vs all-membership period stats)

## [Unreleased] — Script 11 compact ByCluster / ByMouse tables

### Added
- `DominantPeriod_ByCluster.csv`: Cluster, N, Median_h, Median_SD, Mean_h, Mean_SD
- `DominantPeriod_ByMouse.csv`: Cluster, BandName, Median_h, Median_SD, Mean_h, Mean_SD, SignalID, Sex
- Same sheets in `DominantPeriod_Output.xlsx` (`ByCluster`, `ByMouse`)

## [Unreleased] — Script 11 violin only if n≥3

### Changed
- Script 11: KDE violins only when n≥3 points; n<3 shows points + mean/median ticks only

## [Unreleased] — Script 10 panel letters off + LL photoperiod label

### Changed
- All Fig02–07 standalones/composites: A/B/C/E panel lettering disabled (temporary)
- Photoperiod 24 h tick label: `LL` not `L24` (Fig03 co-expression and all `script10_pp_label_` uses)

## [Unreleased] — Script 11 sex-split violin spacing

### Fixed
- `PopulationByCluster_BySex`: F|M violins no longer touch (revert to dx=0.30, violinWidth≤0.22)

## [Unreleased] — Shared cluster colours (Scripts 9–11) + smooth coherence lines

### Added
- `extended_cluster_colour.m` + `pal.cluster` map (Tol Bright): UR_1_3 C01 cyan, UR_3_6 C01 purple, UR_3_6 C02 red

### Changed
- Script 10/9 activity mean lines match coherence twin colour per cluster (not fixed blue / LD red)
- Script 10/9 coherence mean traces: line-only (no markers)
- Script 11 pooled population violins use cluster colours (sex-split panels unchanged: green/yellow)

## [Unreleased] — Script 12 v1.8: clip ≥0, Tukey outliers, dual figures/stats

### Added
- `extended_grouped_x_layout.m` — shared F|M side-by-side spacing (also used in Script 11 population-by-sex)
- Amplitude `Value_Raw` + clipped `Value = max(0, peak − baseline)`
- Tukey 1.5×IQR outlier flags when cell n ≥ 6 (`OutlierExcluded`, `OutlierRule`, `CellN`)
- Primary stats/figures (all data; hollow outliers) + sensitivity (`*_Sens`; outliers omitted)
- Tables: `Sex_Amplitude_Stats_BH_FDR_Primary.csv`, `Sex_Amplitude_Stats_BH_FDR_Sensitivity.csv`, `Sex_Amplitude_Outliers_Tukey.csv`

### Changed
- Box whiskers/geometry clamped at y = 0; non-overlapping F|M boxes via grouped x layout

## [Unreleased] — Script 12 settings fix + amplitude boxplots

### Fixed
- Settings table row mismatch (`AmpMetrics` now single joined string)

### Changed
- Amplitude figures: F|M box-and-whisker side-by-side per photoperiod (replaces mean±SD line plot); stars above each pair

## [Unreleased] — Script 12 pre-LD amplitude (anticipatory)

### Added
- Pre-LD first-peak amplitude (Option A): baseline `[off−W−pre, off−W)`, search `[off−W, off)` — mirrors post-LD
- Figures: `Sex_Amplitude_PreLD_{cluster}` for UR_1_3 C01, UR_3_6 C01, UR_3_6 C02
- Stats: Female vs Male only (no pre-vs-post test); BH-FDR within `ClusterID|Metric` (PreLD and PostLD corrected separately)

## [Unreleased] — Script 12 amplitude ylim floor at 0

### Changed
- Script 12 post-LD amplitude figures: y-axis starts at 0

## [Unreleased] — Script 11 ylim floor at 0

### Changed
- All Script 11 figure y-axes start at 0

## [Unreleased] — Script 11 stable mouse order + pooled population

### Changed
- Per-mouse cluster figures: stable SignalID order (no sort by sex or median); F/M colour only
- Removed `*_ByMedian.*` mouse figures
- Population figures: `PopulationByCluster` = pooled (no sex split); `PopulationByCluster_BySex` = F|M (spaced)

## [Unreleased] — Script 11 violin spacing + ByMedian figures

### Changed
- Population F|M violins: narrower width + larger offset so they no longer touch
- Added per-cluster `*_ByMedian.*` figures (sort by median τ only; sex colouring retained)
- Default per-cluster figures remain sorted by sex, then median τ

## [Unreleased] — Script 11 proper violins + sex colouring

### Changed
- KDE violins taper to pointed tips (support extended ~4× bandwidth; endpoints pinned to zero)
- Sex from SignalID (Script 12 rules): mouse figures coloured F/M; population figure = F|M pair per cluster
- Tables: `Sex` on per-mouse; F/M n + median/mean τ on cluster summary
- Light jitter of raw candidate points on mouse violins; south legend includes Female/Male

## [Unreleased] — Script 11 violin + legend polish

### Changed
- Script 11 violins: smooth symmetric KDE (was jagged/half histogram polygons)
- Removed full-width cohort mean/median guide lines from per-cluster mouse figures
- South-outside legend: magenta dotted = Median, black solid = Mean (all Script 11 figures)

## [Unreleased] — Script 11 separate cluster figs + population panel

### Changed
- Script 11 figures: one file per cluster (`Supp_DominantPeriod_ClusterPeriod_*_C##`)
- Added `Supp_DominantPeriod_PopulationByCluster` — x = UR 1-3 C01 / UR 3-6 C01 / UR 3-6 C02; violin = per-mouse median period spread

## [Unreleased] — Script 11 without WP_TS (mean + median)

### Changed
- Script 11 reads Script 7 `ClusterMembership.RawPeriod_h` only — **no WP_TS load**
- Reports **mean and median** period per mouse and cohort (`MeanPeriod_h` / `MedianPeriod_h`, `MeanTau_Cohort_h` / `MedianTau_Cohort_h`)
- Figure stem: `Supp_DominantPeriod_ClusterPeriod_3Cluster` (candidate-period violins; cohort mean/median guides)

## [Unreleased] — Script 12 amplitude mean±SD, clean panels

### Changed
- Amplitude error bars: mean ± SD (was SEM)
- Removed on-figure footnotes (n F/M, star key); BH-FDR stars retained
- Activity grids: removed per-panel `n =` annotations (n remains in `Sex_N_ByCluster`)

## [Unreleased] — Script 12 amplitude-only (first peak post-LD)

### Changed
- Amplitude metric: first local max after lights-off − pre-off baseline (W ≈ 1× cluster period, cap 3 h; window-max fallback)
- Y-axis: `Amplitude (post-light transition; au)`; figures `Sex_Amplitude_PostLD_*`
- Removed ridge-power panels and day-mean ridge metric; coherence remains out of scope

## [Unreleased] — Script 12 sex amplitude F vs M tests

### Changed
- Amplitude panels: removed “exploratory — not confirmatory” banner
- F vs M two-sided Mann–Whitney U per photoperiod; Cliff’s δ + bootstrap 95% CI; BH-FDR within ClusterID×Metric
- Stars on amplitude figures from Q_BH (* / ** / ***); table `Sex_Amplitude_Stats_BH_FDR.csv`

## [Unreleased] — Script 12 sex-stratified activity + amplitude

### Added
- **Extended Script 12** — Female|Male side-by-side 24h ZT activity grids (L12–L22) for UR_1_3 C01, UR_3_6 C01, UR_3_6 C02
- Amplitude panels (F/M overlaid): mean ridge power + LD-transition spike vs photoperiod; `Sex_Amplitude_PerMouse.csv`
- Entry: `run_extended_script12_sex_stratified_profiles(cohortRoot)`

### Changed
- Removed sex-split phase coherence / WP_TS rebuild (Script 7 `PhaseCoherence_24h` is pooled; coherence omitted from Script 12)

## [Unreleased] — Script 11 dominant UR period summary

### Added
- **Extended Script 11** — read-only dominant validated UR period tables + supplemental 3-panel ridge-period violins (UR_1_3 C01 | UR_3_6 C01 | UR_3_6 C02)
- Entry: `run_extended_script11_dominant_periods(cohortRoot)`; outputs `Script11_DominantPeriod_{Publication|Development}/`
- `DominantPeriod_ClusterSummary.csv` / `DominantPeriod_Output.xlsx` — cohort median τ, IQR across mice, within-mouse IQR, harmonic distance vs P₀ reference
- Paul Tol Bright band colours (`UR_1_3` cyan, `UR_3_6` purple); median (dotted) + mean (solid) overlays
- Script 11 outputs are **not** bundled into Script 10 `Package_Manuscript/` (separate supp deliverable)

## [Unreleased] — Script 10 manuscript amendments

### Changed
- **Script 10** exports comprehensive **supplemental tables** to `Package_Manuscript/Tables/` (`Workbooks/`, `CSV/`, index) and mirrors to `Script10_*/Tables/Supplemental/`; flat CSV aliases kept at `Tables/` root
- Figs 5/7 **main composites = panel A only** (24h ZT coherence grid); transition supp = **LD Pre/Post R** bars (`Supp_Fig05/07_LD_PrePost_R`) under `Standalone/Supp_Transitions/` — not glued into Wide/Tall
- Dropped Script 10 export of `Supp_Fig05/07_B_RidgePower` and `_C_DeltaR` gradient figures (Script 5 FDR tables remain in Package for citation)
- UR_3_6 **C02** (ClusterRank 2, ~4.25–5.0 h) activity + coherence twins as Supp extras (`Supp_Activity_UR36_C02`, `Supp_Coherence_UR36_C02`); mains stay C01
- `cfg.script10.manuscriptClusters` lock: resolve by ClusterRank (+ optional period-window warn); builders accept explicit ClusterID
- Package_Manuscript: main Figs 3–7 unchanged set; Supp under `Figures/Supp/`; Script 5 FDR tables still packaged for text citation
- Figs 5/7 panel A: 24 h ZT phase coherence (Script 7 profiles) instead of peri-transition DL/LD curves
- Activity + coherence grids: figure-level `n = … mice (across photoperiods)` plus per-facet `n = …` (cluster membership can be sparse at early L; union ≠ per panel)
- Dark-phase shade + lights-off line aligned with Script 9 / collaborator style

### Fixed
- `cfg.plotMode = 'publication'` now writes `Script10_ManuscriptFigures_Publication` (mode string normalised in runner + `make_output_dirs_`)
- Fig04/Fig06 (and coherence) `n` mice: count unique `SignalID`/`MouseID` across photoperiods — not `File|SignalID` (File differs per L12–L22 and inflated n above study N=24)
- Fig04/Fig06 annotation clarified: bare `n = 24 mice` was a facet union and misleading when early photoperiods plot fewer grey traces (UR_3_6 C01 membership)

## [Unreleased] — Script 7 activity column match (v2.1)

### Fixed
- `find_table_column` in `extended_script7_profiles_run.m`: sex-prefixed SignalIDs (`F_`/`M_`) now strip before matching Excel headers; fuzzy `contains` no longer uses R2025b scalar broadcasting that bound every mouse to `Time (hr)`. Non-activity columns (time/light/ZT) are excluded. Activity z-scored 24h profiles should again show grey individuals + ultradian structure.

## [Unreleased] — Repo layout cleanup

### Changed
- Legacy monolithic scripts moved to `Legacy/` (`harmonic_subtraction_v6.m`, wavelet v4/v5, `PeriodComparison_v4.m`, `Plot_v5.m`)
- Core docs moved into module folders: `Analysis/METHODS.md`, `Functions/io/DATA_STRUCTURE.md`, `Functions/handoff/HANDOFF_SCHEMA.md`, `Functions/hsub/README.md`
- Root retains project-level docs only (`README.md`, `PROJECT_CONTEXT.md`, `CHANGELOG.md`, `TODO.md`, `CURSOR_MATLAB_WORKFLOW.md`, `LabNotebook.md`)
- `setup_paths()` adds `Legacy/` to the MATLAB path

## [Unreleased] — Extended UR E0–E1

### Added
- Top-level `Extended/` tree with Kent A–G copies in `Legacy/Kent_AG/`
- `setup_extended_paths.m`, `Config/extended_defaults.m`, `EXTENDED_CONTEXT.md`, `METHODS_EXTENDED.md`
- Thin E1 wrappers `Analysis/run_extended_*.m` (A–G via `run()`)

### Notes
- Core Scripts 1–3 frozen; CarryForward ±15% SEL_P360; LL projected dark = L22 ZT22–24
- Core Script 3 residual-CWT vs Extended C are parallel claims

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
- `setup_paths` fprintf notes optional `setup_extended_paths`

### Pending (core-v1.0 tag)
- MATLAB checkpoints Phases 2–6

## [Pre-Core] — Legacy monolithic scripts

- Scripts now in `Legacy/`: `harmonic_subtraction_v6.m`, `behav_wavelet_scalogram_v5.m`, `behav_wavelet_powerspectrum_v4.m`, `PeriodComparison_v4.m`, `Plot_v5.m`
