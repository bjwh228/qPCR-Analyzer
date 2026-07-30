"""
qPCR Import Module
===================
Parses raw qPCR export files from various sources into a normalized
"well-level" long format:

    Sample, Well, Gene, Ct, Group (Group is often absent at this stage)

Supports:
  - Generic long CSV/Excel: Sample, Group, Gene, Ct  (already-clean format)
  - Generic wide CSV/Excel: Sample [,Group], Gene1, Gene2, ...
        or with replicate suffixes: Sample, GAPDH_1, GAPDH_2, IL6_1, IL6_2
  - Bio-Rad CFX Manager exports (Well, Fluor, Target, Sample, Cq, ...)
  - Applied Biosystems / QuantStudio exports (Well Position, Sample Name,
    Target Name, CT, ...) -- these have several metadata header rows before
    the real data table, which this module auto-detects and skips.

Missing-value markers ("Undetermined", "N/A", "No Ct", blank, etc.) are
normalized to NaN here; what to DO about them is a separate, user-chosen
decision handled in qpcr_qc.py.
"""

import pandas as pd
import numpy as np
import io
import re

MISSING_MARKERS = {
    "undetermined", "undet", "n/a", "na", "no ct", "nocall", "unknown",
    "none", "nan", "", "-", "no amp", "no amplification", "fail", "failed",
}

# column name -> list of aliases seen across instruments/software versions
COLUMN_ALIASES = {
    "sample": ["sample", "sample name", "sample id", "sampleid"],
    "gene":   ["gene", "target", "target name", "assay", "assay name", "fluor"],
    "ct":     ["ct", "cq", "cq value", "ct value", "c(t)", "c t", "c\u0442"],  # last handles mangled encoding
    "well":   ["well", "well position", "well pos"],
    "group":  ["group", "condition", "treatment", "biological group"],
}


def _normalize_col(name):
    return re.sub(r"\s+", " ", str(name).strip().lower())


def _match_alias(colname):
    norm = _normalize_col(colname)
    for standard, aliases in COLUMN_ALIASES.items():
        if norm in aliases:
            return standard
    return None


def _is_missing(val):
    if pd.isna(val):
        return True
    return str(val).strip().lower() in MISSING_MARKERS


def _to_ct_float(val):
    if _is_missing(val):
        return np.nan
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def _find_quantstudio_header_row(lines, max_scan=60):
    """
    QuantStudio exports reliably have a data header row starting with 'Well'
    (e.g. 'Well,Sample Name,Target Name,Task,CT,...' or 'Well Position,...').
    If the file's first line already starts with 'Well', no metadata rows
    precede it. Otherwise, scan forward for the first line that does and
    treat that as the header, skipping everything above it (QuantStudio
    prints several '* Block Type: ...' style metadata lines before the table).
    """
    for i, line in enumerate(lines[:max_scan]):
        if line.strip().lower().startswith("well"):
            return i
    return None


def _find_header_row(lines, max_scan=60):
    """Scan raw text lines for the row that looks like the real data header
    (matches >=2 known column aliases), skipping instrument metadata rows."""
    quantstudio_idx = _find_quantstudio_header_row(lines, max_scan=max_scan)
    if quantstudio_idx is not None:
        return quantstudio_idx

    best_idx, best_score = None, 0
    for i, line in enumerate(lines[:max_scan]):
        cells = [c.strip() for c in re.split(r",|\t", line)]
        score = sum(1 for c in cells if _match_alias(c) is not None)
        if score > best_score:
            best_score, best_idx = score, i
    return best_idx if best_score >= 2 else 0


def _rename_with_aliases(df):
    """
    Map raw column names to standard names (sample/gene/ct/well/group).
    If more than one raw column maps to the same standard name (e.g. a
    QuantStudio export with both 'Well' and 'Well Position'), keep only the
    first match and drop the rest -- otherwise pandas silently creates
    duplicate-named columns that break downstream selection.
    """
    rename_map = {}
    used_standards = set()
    cols_to_drop = []
    for col in df.columns:
        std = _match_alias(col)
        if std:
            if std in used_standards:
                cols_to_drop.append(col)
            else:
                rename_map[col] = std
                used_standards.add(std)
    df = df.drop(columns=cols_to_drop)
    return df.rename(columns=rename_map)


def _melt_wide_format(df):
    """
    Handles wide format where each gene (optionally with replicate suffix
    like GAPDH_1, GAPDH_2, or GAPDH.1 from pandas auto-dedup) is its own
    column, e.g.:
        Sample, Group, GAPDH_1, GAPDH_2, IL6_1, IL6_2
    """
    id_cols = [c for c in ["sample", "group", "well"] if c in df.columns]
    value_cols = [c for c in df.columns if c not in id_cols]

    melted = df.melt(id_vars=id_cols, value_vars=value_cols,
                      var_name="gene_raw", value_name="ct")

    # strip trailing replicate markers like "_1", ".1", " rep2", "-2"
    melted["gene"] = melted["gene_raw"].str.replace(
        r"([_.\-\s]?(rep)?\d+)$", "", regex=True, case=False
    ).str.strip()
    melted = melted.drop(columns=["gene_raw"])
    return melted


def load_qpcr_file(file_obj, filename=None):
    """
    Main entry point. Accepts a file path (str), a file-like object, or bytes.
    Returns (long_df, format_detected) where long_df has columns:
        Sample, Gene, Ct, [Well], [Group]
    Note: rows may contain duplicate (Sample, Gene) pairs -- these are
    technical replicates and should be passed to qpcr_qc.collapse_replicates()
    before running ddCt/statistics.
    """
    # ---- read raw bytes/text ----
    is_excel = False
    if filename:
        is_excel = filename.lower().endswith((".xlsx", ".xls"))
    elif hasattr(file_obj, "name"):
        is_excel = file_obj.name.lower().endswith((".xlsx", ".xls"))

    if is_excel:
        raw_df = pd.read_excel(file_obj, header=None)
        lines = [",".join(str(x) for x in row) for row in raw_df.values.tolist()]
    else:
        if isinstance(file_obj, bytes):
            text = file_obj.decode("utf-8", errors="replace")
        elif hasattr(file_obj, "read"):
            content = file_obj.read()
            text = content.decode("utf-8", errors="replace") if isinstance(content, bytes) else content
        else:
            with open(file_obj, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        lines = text.splitlines()

    header_idx = _find_header_row(lines)

    if is_excel:
        df = pd.read_excel(file_obj, header=header_idx)
    else:
        # detect delimiter from the header line itself
        delim = "\t" if lines[header_idx].count("\t") > lines[header_idx].count(",") else ","
        # skip_blank_lines=False keeps row indices aligned with _find_header_row,
        # which scanned raw (unfiltered) lines -- otherwise blank metadata lines
        # shift the header row index pandas actually uses.
        df = pd.read_csv(io.StringIO(text), header=header_idx, sep=delim,
                          engine="python", skip_blank_lines=False)
        df = df.dropna(how="all")

    df.columns = [str(c) for c in df.columns]
    df = _rename_with_aliases(df)

    format_detected = "unknown"

    if {"sample", "gene", "ct"}.issubset(df.columns):
        format_detected = "long"
        keep = [c for c in ["sample", "group", "gene", "ct", "well"] if c in df.columns]
        long_df = df[keep].copy()
    elif "sample" in df.columns:
        format_detected = "wide"
        long_df = _melt_wide_format(df)
    else:
        raise ValueError(
            "Could not detect a 'Sample' column (or alias like 'Sample Name'). "
            f"Columns found: {list(df.columns)}. If this is a raw instrument export "
            "we haven't seen before, share a few example rows and the parser can be extended."
        )

    long_df["ct"] = long_df["ct"].apply(_to_ct_float)
    long_df = long_df.rename(columns={
        "sample": "Sample", "gene": "Gene", "ct": "Ct", "well": "Well", "group": "Group"
    })
    long_df["Sample"] = long_df["Sample"].astype(str).str.strip()
    long_df["Gene"] = long_df["Gene"].astype(str).str.strip()

    return long_df.reset_index(drop=True), format_detected


if __name__ == "__main__":
    # smoke test: generic long CSV
    csv_text = "Sample,Group,Gene,Ct\nS1,Control,GAPDH,18.2\nS1,Control,IL6,24.1\n"
    df, fmt = load_qpcr_file(io.StringIO(csv_text).read().encode(), filename="test.csv")
    print("Long format test:\n", df, "\nDetected:", fmt)

    # smoke test: wide format with replicate suffixes
    csv_text2 = "Sample,Group,GAPDH_1,GAPDH_2,IL6_1,IL6_2\nS1,Control,18.1,18.3,24.0,24.2\n"
    df2, fmt2 = load_qpcr_file(io.StringIO(csv_text2).read().encode(), filename="test2.csv")
    print("\nWide format test:\n", df2, "\nDetected:", fmt2)

    # smoke test: mock QuantStudio-style export with metadata header rows + Undetermined
    csv_text3 = (
        "* Block Type: 96\n* Instrument: QuantStudio\n\n"
        "Well Position,Sample Name,Target Name,CT\n"
        "A1,S1,GAPDH,18.4\n"
        "A2,S1,IL6,Undetermined\n"
    )
    df3, fmt3 = load_qpcr_file(io.StringIO(csv_text3).read().encode(), filename="test3.csv")
    print("\nQuantStudio-style test:\n", df3, "\nDetected:", fmt3)