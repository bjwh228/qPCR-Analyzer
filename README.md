# qPCR Analyzer

A qPCR analysis tool that takes raw Ct values, computes ddCt fold-change,
runs the appropriate statistical test automatically, and generates a publication-ready,
fully customizable figure in Python or R.

## Files

| File | Purpose |
|---|---|
| `qpcr_import.py` | Parses multiple file formats into a normalized well-level table |
| `qpcr_qc.py` | Missing-value handling + technical replicate QC |
| `qpcr_core.py` | ddCt calculation + statistics |
| `qpcr_plot.py` | Customizable plotting functions |
| `app.py` | Streamlit web app |
| `qpcr_analysis.R` | Standalone R/ggplot2 CLI script |
| `sample_data.csv` | Example input |
| `example_instrument_export.csv` | Emulates a QuantStudio-style export with common QC issues |

## Supported input formats

The importer (`qpcr_import.py`) auto-detects and handles:

- **Clean long CSV/Excel**: `Sample, Group, Gene, Ct` (one row per measurement)
- **Wide CSV/Excel**: `Sample, Group, GAPDH_1, GAPDH_2, IL6_1, IL6_2, ...` (auto-melted, replicate suffixes stripped to recover the gene name)
- **Instrument exports** (Bio-Rad CFX Manager, Applied Biosystems/QuantStudio): scans for the header row automatically by matching known column-name aliases (`Sample Name`, `Target Name`, `CT`/`Cq`, `Well Position`, etc.)
- Missing/failed calls are recognized in any of these forms: `Undetermined`, `N/A`, `No Ct`, blank, `Unknown`, etc., and normalized to a proper missing value.

More upload formats will be added in the future.

If a file has no `Group` column, an editable table appears where you assign each sample to a group before proceeding.

## Handling duplicates (technical replicates) and missing values

Common qPCR protocols have technical replicate wells and occasional undertermined expression levels. `qpcr_qc.py` handles both, with user-facing choices in the app:

**Missing/undetermined Ct values**:
- *Exclude just that well* (default) — keep other valid replicates for that Sample+Gene
- *Impute a fixed Ct* — replace with a chosen value (e.g. 40) representing "no amplification"
- *Exclude the whole Sample+Gene* if any replicate is missing (most conservative)

**Technical replicate agreement** — standard qPCR QC practice is to check
that replicate wells agree within a tolerance before averaging them. You set:
- **Cutoff**, in either fold-change (e.g. flag if replicates differ by ≥2 fold, which equals a 1 Ct spread) or raw Ct units (e.g. flag if they span ≥0.5 Ct)
- **Handling Disagreeing Values**:
  - *Flag* — averages everything, but marks it for review in the QC Report tab (default)
  - *Exclude* — drops that Sample+Gene data point entirely
  - *Drop* — with 3+ replicates, removes the single most-deviant value and re-averages the rest (falls back to exclude with only 2 replicates)

A log is generated from **QC Report** tab to aid in reproducibility.

## Input format

This format works best for analysis:

```csv
Sample,Group,Gene,Ct
Control_1,Control,GAPDH,18.20
Control_1,Control,IL6,18.06
Treatment_1,Treatment,GAPDH,18.15
Treatment_1,Treatment,IL6,21.30
```

## Notes on the statistics

- Tests are run on **dCt values, not fold-change** — fold-change (2^-ddCt) is
  log-distributed, so running a t-test/ANOVA directly on it violates the
  normality assumption. This is a common mistake in qPCR analysis.

- Tukey HSD is used as the standard post-hoc test after a significant
  one-way ANOVA; this can be swapped for Dunnett's test
