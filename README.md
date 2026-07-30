# qPCR Analyzer

A tool that takes raw Ct values, computes ddCt fold-change (Livak method),
runs the appropriate statistical test automatically (Welch's t-test for 2
groups, one-way ANOVA + Tukey HSD for 3+), and generates a publication-ready,
fully customizable figure — in Python (matplotlib/seaborn, via a Streamlit
web app) or R (ggplot2, via command line).

## Files

| File | Purpose |
|---|---|
| `qpcr_import.py` | Parses multiple file formats into a normalized well-level table |
| `qpcr_qc.py` | Missing-value handling + technical replicate QC (flag/exclude/drop-outlier) |
| `qpcr_core.py` | ddCt calculation + statistics engine (no UI) |
| `qpcr_plot.py` | Customizable matplotlib/seaborn plotting functions |
| `app.py` | Streamlit web app — the actual product, upload-and-go |
| `qpcr_analysis.R` | Standalone R/ggplot2 CLI script, same method (no import/QC layer yet) |
| `sample_data.csv` | Clean example input, already-averaged Ct values |
| `example_instrument_export.csv` | Messy example: mimics a QuantStudio-style export with metadata header rows, an `Undetermined` value, and disagreeing technical replicates — good for testing the QC pipeline |

## Supported input formats

The importer (`qpcr_import.py`) auto-detects and handles:

- **Clean long CSV/Excel**: `Sample, Group, Gene, Ct` (one row per measurement)
- **Wide CSV/Excel**: `Sample, Group, GAPDH_1, GAPDH_2, IL6_1, IL6_2, ...` (auto-melted, replicate suffixes stripped to recover the gene name)
- **Instrument exports** (Bio-Rad CFX Manager, Applied Biosystems/QuantStudio style): these have several metadata rows before the real header — the importer scans for the header row automatically by matching known column-name aliases (`Sample Name`, `Target Name`, `CT`/`Cq`, `Well Position`, etc.)
- Missing/failed calls are recognized in any of these forms: `Undetermined`, `N/A`, `No Ct`, blank, `Unknown`, etc., and normalized to a proper missing value.

If your instrument's export doesn't parse cleanly, share a few example rows and the alias/header-detection logic in `qpcr_import.py` can be extended — this is the piece most likely to need a small tweak per lab, since every instrument's software version formats slightly differently.

If a file has no `Group` column (common for raw instrument exports — the machine doesn't know your experimental design), the Streamlit app shows an editable table where you assign each sample to a group before proceeding.

## Handling duplicates (technical replicates) and missing values

Real qPCR data always has technical replicate wells and occasional failed
calls. `qpcr_qc.py` handles both, with user-facing choices in the app:

**Missing/undetermined Ct values** — three strategies:
- *Exclude just that well* (default) — keep other valid replicates for that Sample+Gene
- *Impute a fixed Ct* — replace with a chosen value (e.g. 40) representing "no amplification"
- *Exclude the whole Sample+Gene* if any replicate is missing (most conservative)

**Technical replicate agreement** — standard qPCR QC practice is to check
that replicate wells agree within a tolerance before averaging them. You set:
- **Cutoff**, in either fold-change (e.g. flag if replicates differ by ≥2-fold, which equals a 1-Ct spread) or raw Ct units (e.g. flag if they span ≥0.5 Ct — a common core-facility threshold)
- **What to do when they disagree**:
  - *Flag but keep* — averages everything, but marks it for your review in the QC Report tab (default; doesn't silently discard data)
  - *Exclude* — drops that Sample+Gene data point entirely
  - *Drop the outlier* — with 3+ replicates, removes the single most-deviant value and re-averages the rest (falls back to exclude with only 2 replicates, since there's no way to tell which of two disagreeing values is the "bad" one)

Every decision is logged and downloadable from the app's **QC Report** tab, so nothing is silently dropped without a paper trail — important for methods sections and reproducibility.

## Input format (clean data / R script)

If your data is already clean and averaged, this format still works directly:

```csv
Sample,Group,Gene,Ct
Control_1,Control,GAPDH,18.20
Control_1,Control,IL6,18.06
Treatment_1,Treatment,GAPDH,18.15
Treatment_1,Treatment,IL6,21.30
```

The R script (`qpcr_analysis.R`) currently expects this clean format directly — the import/QC layer is Python-only for now.

## Run it locally (Python / Streamlit)

```bash
pip install streamlit pandas numpy scipy statsmodels matplotlib seaborn
streamlit run app.py
```

This opens a browser UI where users upload their CSV, pick the reference
gene / control group, tweak colors/fonts/figure size, and download the
figure (PNG/PDF/SVG) plus stats tables (CSV).

## Run it in R

```bash
Rscript qpcr_analysis.R sample_data.csv --reference GAPDH --control Control --output figure.pdf
```

## Getting this to revenue quickly

**Fastest path to a live, payable product:**

1. **Deploy the Streamlit app for free** on [Streamlit Community Cloud](https://streamlit.io/cloud) or Render's free tier — takes ~10 minutes, gives you a public URL immediately.
2. **Add a simple usage gate.** Two low-effort options:
   - *Freemium*: free for e.g. 1 gene / 2 groups, paywall multi-gene or 3+ group ANOVA analyses.
   - *Pay-per-report*: Stripe Payment Links (no code) gating a "download high-res figure" button — a $2–5 one-time charge per analysis is a very low-friction ask for a grad student who just needs one clean figure for a paper/thesis.
3. **Distribution (all free, all fast):**
   - Post in r/labrats, r/bioinformatics, r/AskAcademia with a before/after screenshot (messy Excel chart vs. your clean, annotated bar chart).
   - Bluesky/X science community, university core-facility Slack/Discord channels.
   - A short before/after post ("stop making qPCR bar charts in Excel — do this in 30 seconds") tends to do well because the pain point is universal and the visual is self-explanatory.
4. **Upsell path once you have a few users:** batch analysis (many genes/plates at once), an "export a full stats + methods paragraph for your paper" feature (huge time-saver, easy to charge more for), and journal-specific style presets (Nature, Cell, PLOS figure formatting conventions).

**Why this is more likely to convert than a bigger tool:** the buyer is an
individual researcher with a credit card and an urgent, well-defined task
(I need this figure by Friday), not an institution that needs a purchase
order. Keep the free tier genuinely useful so the "aha, this saved me an
hour" moment happens before you ask for money.

## Notes on the statistics

- Tests are run on **dCt values, not fold-change** — fold-change (2^-ddCt) is
  log-distributed, so running a t-test/ANOVA directly on it violates the
  normality assumption. This is a common mistake in home-grown Excel
  templates and is worth mentioning in your marketing (it signals rigor to
  a scientific audience).
- Tukey HSD is used as the standard post-hoc test after a significant
  one-way ANOVA; this can be swapped for Dunnett's test (compare everything
  to control only, not all pairs) if that fits your users' needs better —
  worth A/B testing as a feature.
