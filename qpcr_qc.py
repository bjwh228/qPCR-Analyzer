"""
qPCR QC Module
==============
Two jobs, applied in order to well-level data (which may contain duplicate
technical-replicate rows per Sample+Gene):

  1. handle_missing_values() -- decide what to do with Undetermined/NaN Ct
     values, per user-chosen strategy.
  2. collapse_replicates()   -- average technical replicates down to one Ct
     per (Sample, Gene), flagging or excluding replicate sets that disagree
     by more than a user-set cutoff (standard qPCR QC practice).
"""

import pandas as pd
import numpy as np


# ---------------------------------------------------------------- missing values

def handle_missing_values(df, ct_col="Ct", strategy="exclude_well", impute_ct=40.0):
    """
    Parameters
    ----------
    strategy : one of
        "exclude_well"      -- drop just the individual missing well/replicate;
                                remaining replicates for that Sample+Gene are
                                still used (recommended default)
        "impute_fixed_ct"   -- replace missing Ct with a fixed value (e.g. 40,
                                representing the assay's no-amplification cutoff)
        "exclude_sample_gene" -- if ANY replicate for a Sample+Gene is missing,
                                drop the entire Sample+Gene data point (most
                                conservative)
    Returns
    -------
    (clean_df, log_df) -- log_df records every row that was dropped/imputed and why
    """
    df = df.copy()
    missing_mask = df[ct_col].isna()
    log_rows = []

    if strategy == "exclude_well":
        for _, row in df[missing_mask].iterrows():
            log_rows.append({**row.to_dict(), "action": "excluded (missing well)"})
        clean = df[~missing_mask].copy()

    elif strategy == "impute_fixed_ct":
        for _, row in df[missing_mask].iterrows():
            log_rows.append({**row.to_dict(), "action": f"imputed Ct={impute_ct} (missing/undetermined)"})
        clean = df.copy()
        clean.loc[missing_mask, ct_col] = impute_ct

    elif strategy == "exclude_sample_gene":
        bad_pairs = df[missing_mask][["Sample", "Gene"]].drop_duplicates()
        merged = df.merge(bad_pairs, on=["Sample", "Gene"], how="left", indicator=True)
        to_drop = merged["_merge"] == "both"
        for _, row in df[to_drop.values].iterrows():
            log_rows.append({**row.to_dict(), "action": "excluded (Sample+Gene has a missing replicate)"})
        clean = df[~to_drop.values].copy()

    else:
        raise ValueError(f"Unknown missing-value strategy: {strategy}")

    log_df = pd.DataFrame(log_rows) if log_rows else pd.DataFrame(
        columns=list(df.columns) + ["action"])
    return clean.reset_index(drop=True), log_df


# ---------------------------------------------------------------- replicate collapsing

def collapse_replicates(df, ct_col="Ct", cutoff=1.0, cutoff_type="fold_change",
                         on_disagreement="flag"):
    """
    Collapse technical replicates (duplicate Sample+Gene rows) into one mean
    Ct per Sample+Gene, flagging or excluding replicate sets that disagree
    by more than `cutoff`.

    Parameters
    ----------
    cutoff_type : "fold_change" (default) or "ct"
        "fold_change": cutoff is expressed as fold-difference, e.g. cutoff=1.0
            means flag if replicates differ by 2^1.0 = 2-fold or more.
            (Ct difference of 1.0 = exactly 2-fold; this is the intuitive
            unit for biologists, hence the default.)
        "ct": cutoff is a raw Ct-value difference, e.g. cutoff=0.5 flags
            replicate sets whose Ct values span more than 0.5 cycles
            (a common core-facility QC threshold).
    on_disagreement : "flag" | "exclude" | "drop_outlier"
        "flag"        -- keep the mean of all replicates but mark as flagged
                         so it's visible in the QC report (recommended default,
                         doesn't silently throw away data)
        "exclude"     -- drop the entire Sample+Gene data point
        "drop_outlier"-- with 3+ replicates, drop the single value farthest
                         from the median and re-average the rest; with only
                         2 replicates that disagree, falls back to "exclude"
                         (no way to tell which of 2 values is the outlier)

    Returns
    -------
    (collapsed_df, qc_report_df)
      collapsed_df : Sample, Group (if present), Gene, Ct, n_replicates,
                     replicate_ct_range, replicate_fold_range, flagged
      qc_report_df : one row per Sample+Gene showing all replicate values,
                     range, and what action was taken -- for full transparency
    """
    if cutoff_type not in ("fold_change", "ct"):
        raise ValueError("cutoff_type must be 'fold_change' or 'ct'")

    group_cols = ["Sample", "Gene"]
    extra_cols = [c for c in ["Group"] if c in df.columns]

    rows = []
    qc_rows = []

    for (sample, gene), sub in df.groupby(group_cols):
        vals = sub[ct_col].dropna().values
        n = len(vals)
        group_val = sub[extra_cols[0]].iloc[0] if extra_cols else None

        if n == 0:
            continue

        if n == 1:
            rows.append({"Sample": sample, "Gene": gene, "Group": group_val,
                         "Ct": vals[0], "n_replicates": 1,
                         "replicate_ct_range": 0.0, "replicate_fold_range": 1.0,
                         "flagged": False})
            qc_rows.append({"Sample": sample, "Gene": gene, "replicate_values": list(vals),
                             "n_replicates": 1, "ct_range": 0.0, "fold_range": 1.0,
                             "flagged": False, "action": "single replicate, no QC check possible"})
            continue

        ct_range = vals.max() - vals.min()
        fold_range = 2 ** ct_range
        is_bad = (fold_range >= (2 ** cutoff)) if cutoff_type == "fold_change" else (ct_range >= cutoff)

        action = "mean of all replicates (within cutoff)"
        final_vals = vals
        flagged = False

        if is_bad:
            flagged = True
            if on_disagreement == "flag":
                action = f"FLAGGED but kept: replicates span {ct_range:.2f} Ct ({fold_range:.2f}-fold)"
            elif on_disagreement == "exclude":
                qc_rows.append({"Sample": sample, "Gene": gene, "replicate_values": list(vals),
                                 "n_replicates": n, "ct_range": ct_range, "fold_range": fold_range,
                                 "flagged": True, "action": "EXCLUDED: replicates disagree beyond cutoff"})
                continue
            elif on_disagreement == "drop_outlier":
                if n >= 3:
                    median = np.median(vals)
                    idx_outlier = np.argmax(np.abs(vals - median))
                    final_vals = np.delete(vals, idx_outlier)
                    new_range = final_vals.max() - final_vals.min() if len(final_vals) > 1 else 0.0
                    action = (f"dropped outlier replicate ({vals[idx_outlier]:.2f}), "
                              f"averaged remaining {len(final_vals)} (new range {new_range:.2f} Ct)")
                    flagged = new_range >= cutoff if cutoff_type == "ct" else (2 ** new_range) >= (2 ** cutoff)
                else:
                    qc_rows.append({"Sample": sample, "Gene": gene, "replicate_values": list(vals),
                                     "n_replicates": n, "ct_range": ct_range, "fold_range": fold_range,
                                     "flagged": True,
                                     "action": "EXCLUDED: only 2 replicates disagree, can't identify outlier"})
                    continue
            else:
                raise ValueError(f"Unknown on_disagreement option: {on_disagreement}")

        mean_ct = final_vals.mean()
        rows.append({"Sample": sample, "Gene": gene, "Group": group_val,
                     "Ct": mean_ct, "n_replicates": len(final_vals),
                     "replicate_ct_range": ct_range, "replicate_fold_range": fold_range,
                     "flagged": flagged})
        qc_rows.append({"Sample": sample, "Gene": gene, "replicate_values": list(vals),
                         "n_replicates": n, "ct_range": ct_range, "fold_range": fold_range,
                         "flagged": flagged, "action": action})

    collapsed = pd.DataFrame(rows)
    qc_report = pd.DataFrame(qc_rows)

    if extra_cols and "Group" in collapsed.columns and collapsed["Group"].isna().all():
        collapsed = collapsed.drop(columns=["Group"])

    return collapsed.reset_index(drop=True), qc_report.reset_index(drop=True)


if __name__ == "__main__":
    # smoke test
    df = pd.DataFrame([
        {"Sample": "S1", "Group": "Control", "Gene": "GAPDH", "Ct": 18.1},
        {"Sample": "S1", "Group": "Control", "Gene": "GAPDH", "Ct": 18.3},
        {"Sample": "S1", "Group": "Control", "Gene": "IL6", "Ct": 24.0},
        {"Sample": "S1", "Group": "Control", "Gene": "IL6", "Ct": 24.2},
        # bad replicate pair: 2-Ct spread (4-fold) on IL6 for S2
        {"Sample": "S2", "Group": "Control", "Gene": "IL6", "Ct": 22.0},
        {"Sample": "S2", "Group": "Control", "Gene": "IL6", "Ct": 24.0},
        # 3 replicates, one outlier
        {"Sample": "S3", "Group": "Treatment", "Gene": "IL6", "Ct": 20.0},
        {"Sample": "S3", "Group": "Treatment", "Gene": "IL6", "Ct": 20.1},
        {"Sample": "S3", "Group": "Treatment", "Gene": "IL6", "Ct": 22.5},
        # missing well
        {"Sample": "S4", "Group": "Treatment", "Gene": "IL6", "Ct": np.nan},
        {"Sample": "S4", "Group": "Treatment", "Gene": "IL6", "Ct": 21.0},
    ])

    print("=== missing value handling (exclude_well) ===")
    clean, log = handle_missing_values(df, strategy="exclude_well")
    print(log[["Sample", "Gene", "action"]])

    print("\n=== collapse_replicates, on_disagreement='flag' ===")
    collapsed, qc = collapse_replicates(clean, cutoff=1.0, cutoff_type="fold_change", on_disagreement="flag")
    print(collapsed)
    print("\nQC report:")
    print(qc[["Sample", "Gene", "n_replicates", "fold_range", "flagged", "action"]])

    print("\n=== collapse_replicates, on_disagreement='drop_outlier' ===")
    collapsed2, qc2 = collapse_replicates(clean, cutoff=1.0, cutoff_type="fold_change", on_disagreement="drop_outlier")
    print(collapsed2)
