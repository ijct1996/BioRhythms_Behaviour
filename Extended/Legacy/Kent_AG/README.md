# Kent A–G Final Versions (archived copies)

Immutable copies of collaborator Kent scripts from:

`…\Kent\05. MATLAB\03. Mice_behav\Codes\Final Versions\`

Do not edit these files in place. Modularisation is planned for Extended E2 via `Extended/Functions/`. Entry points are thin wrappers in `Extended/Analysis/`.

## Mapping A–G → `run_extended_*`

| Kent letter | Legacy file in this folder | Wrapper |
|-------------|----------------------------|---------|
| **A** | `A_Harmonic_subtraction_v12.m` | `run_extended_hsub` |
| **B** | `B_Behav_wavelet_v12.m` | `run_extended_wavelet_ridge` |
| **C** | `C_Validate_RawVsSelectiveHSub_Periods_v1.m` | `run_extended_period_gate` |
| **D** | `D_Ultradian_RidgePhase_Resync_v4_FDR_LLProjected.m` | `run_extended_ridge_resync` |
| **D aux** | `D_RidgePeriodGraphing.m` | (companion; called/noted from D wrapper) |
| **D.5** | `D5_Plotting_Ridge_Periods.m` *(renamed from extensionless `D.5_Plotting_Ridge_Periods`)* | (companion; noted from D wrapper) |
| **E** | `E_AcrossPhotoperiod_Analysis_v8_1_ValidatedRawUR_FDR_SexStable_PCAFix.m` | `run_extended_across_photoperiod` |
| **F** | `F_Plot_ValidatedUR_PhaseEventHistograms_v1.m` | `run_extended_phase_events` |
| **G** | `G_Plot_SelectedValidatedUR_PublicationProfiles_v1.m` | `run_extended_publication_profiles` |

## Typical order

```
A → B → C → D (+ D aux / D.5) → E → F → G
```

## Relation to Core v1

Core Scripts 1–3 (`Analysis/run_harmonic_subtraction`, `run_behav_wavelet`, `run_across_photoperiod`) are **frozen** and scientifically parallel, not upstream of this tree. See `Extended/METHODS_EXTENDED.md` for Core Script 3 residual-CWT vs Extended C (CarryForward) as **parallel claims**.
