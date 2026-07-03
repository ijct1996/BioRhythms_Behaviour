# Phase 6 — Full cohort regression checklist

Run in MATLAB R2025b. **Pause OneDrive sync** before batch.

## Setup

```matlab
cd('C:\Users\User\Dev\Cursor\Research\Chronobiology\BioRhythms_Behaviour')
setup_paths()
copyfile('user_paths.example.m','user_paths.m')  % if needed
```

Edit `user_paths.m`:
```matlab
up.outputRoot = 'C:\ResearchData\Results\Core';
```

## Input root

```
C:\Users\User\Dev\Cursor\Research\Chronobiology\Data\BioRhythms Behaviour Data\Pipeline_Input\
```

## Script 1 + 2 (per cohort)

| Cohort | Files | HSub output | Wavelet output |
|--------|-------|-------------|----------------|
| C57_LP | `C57\L*.xlsx` (7) | `...\C57_LP\01_HarmonicSubtraction\` | `...\C57_LP\02_Wavelet\` |
| NR2B_LP | `NR2B_LP\L*.xlsx` (7) | `...\NR2B_LP\01_HarmonicSubtraction\` | `...\NR2B_LP\02_Wavelet\` |
| NR2B_LD_DD | 2 files | `...\NR2B_LD_DD\01_HarmonicSubtraction\` | `...\NR2B_LD_DD\02_Wavelet\` |

Handoff accumulates under each cohort's `Handoff/` (sibling of `02_Wavelet`).

## Script 3 (once per cohort)

Point at `{cohort}/Handoff/`; tag C57_LP / NR2B_LP / NR2B_LD_DD.

## Pass criteria

- [ ] All Anchor_Report sheets readable; majority AnchorOK at L12
- [ ] Recommendations default SEL_P360 when AnchorOK
- [ ] Jet scalograms for every mouse
- [ ] 7 (or 2) `CoreSummary__*.mat` per cohort
- [ ] Co-expression plots span photoperiod axis
