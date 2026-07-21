# Legacy scripts (pre–Core v1)

These monolithic scripts live in this folder and are superseded by the three-script Core v1 pipeline in `Analysis/`. `setup_paths()` adds `Legacy/` to the MATLAB path for reference runs.

| Legacy file | Replaced by |
|-------------|-------------|
| `Legacy/harmonic_subtraction_v6.m` | `Analysis/run_harmonic_subtraction.m` |
| `Legacy/behav_wavelet_scalogram_v5.m` | `Analysis/run_behav_wavelet.m` (scalogram arm) |
| `Legacy/behav_wavelet_powerspectrum_v4.m` | `Analysis/run_behav_wavelet.m` (period + handoff arm) |
| `Legacy/PeriodComparison_v4.m` | `Analysis/run_across_photoperiod.m` |
| `Legacy/Plot_v5.m` | Plotting in `Functions/plot/` + Script 3 |

Kept for reference until Core v1 regression passes (Phase 6). Do not extend.
