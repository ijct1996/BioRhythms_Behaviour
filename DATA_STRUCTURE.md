# Data Structure: {{PROJECT_NAME}}

> Document every dataset used in this project. **Do not commit raw identifiable patient data to Git.**

## Storage locations

| Dataset | Path (local) | Git | Notes |
|---------|--------------|-----|-------|
| Raw | `Data/raw/` | No (gitignored) | |
| Processed | `Data/processed/` | No | |
| Metadata | `Data/metadata/` | No | |

## Primary dataset: *(name)*

**File**: `Data/raw/...`
**Format**: CSV / MAT / XLSX
**Rows**: samples or observations
**Last updated**:

### Columns

| # | Name | Type | Unit | Missing values | Acceptable range | Description |
|---|------|------|------|----------------|------------------|-------------|
| 1 | patient_id | string | — | none expected | — | Unique patient identifier (anonymised) |
| 2 | sampling_time | datetime | — | | | Time of sample collection |
| 3 | sex | categorical | — | | M, F | Biological sex |
| 4 | age | numeric | years | | 18–100 | Age at sampling |
| | | | | | | |

### Group definitions

| Group label | Definition | N (expected) |
|-------------|------------|--------------|
| Control | | |
| Case | | |

### Coding notes

<!-- e.g. 999 = missing, -1 = below detection limit -->

## Derived variables

| Variable | Formula / derivation | Used in |
|----------|---------------------|---------|
| | | |

## Data provenance

- **Source**:
- **Collection dates**:
- **Exclusions applied before import**:
