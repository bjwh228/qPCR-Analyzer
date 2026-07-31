"""
qPCR Analysis Core Engine
=========================
Computes ddCt fold-change values from raw Ct data and runs appropriate
statistical tests (t-test for 2 groups, one-way ANOVA + Tukey HSD for 3+).

Expected input CSV format (long format):
    Sample, Group, Gene, Ct
    S1, Control, GAPDH, 18.2
    S1, Control, TargetGeneA, 24.1
    S2, Control, GAPDH, 18.4
    S2, Control, TargetGeneA, 24.3
    S3, Treatment, GAPDH, 18.1
    S3, Treatment, TargetGeneA, 21.0
    ...

- One row per (Sample, Gene) Ct measurement.
- 'Gene' must include your reference/housekeeping gene(s) and target gene(s).
- Each Sample must belong to exactly one Group.
- If you have technical replicate wells per sample, average them to a single
  Ct per (Sample, Gene) BEFORE feeding this tool (biological replicates are
  the unit of statistical analysis, not technical replicates).
"""

import pandas as pd
import numpy as np
from scipy import stats as sp_stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import dunnett

def compute_dct(df, reference_gene, target_genes=None,
                 sample_col="Sample", group_col="Group", gene_col="Gene", ct_col="Ct"):
    """
    Compute delta-Ct only (Ct_target - Ct_reference), with no control-group
    normalization or fold-change step. Useful when the user hasn't defined
    experimental groups yet and just wants QC'd, normalized Ct values.

    Returns a DataFrame with: Sample, [Group if present], Gene, Ct, Ct_ref, dCt
    """
    df = df.copy()
    if target_genes is None:
        target_genes = [g for g in df[gene_col].unique() if g != reference_gene]

    ref = df[df[gene_col] == reference_gene][[sample_col, ct_col]].rename(
        columns={ct_col: "Ct_ref"})
    if ref[sample_col].duplicated().any():
        raise ValueError(
            f"Multiple '{reference_gene}' Ct values found for the same sample. "
            "Average technical replicates before running this tool."
        )

    results = []
    for gene in target_genes:
        sub = df[df[gene_col] == gene].merge(ref, on=sample_col, how="left")
        if sub["Ct_ref"].isna().any():
            missing = sub[sub["Ct_ref"].isna()][sample_col].tolist()
            raise ValueError(f"No reference-gene Ct found for sample(s): {missing}")
        sub["dCt"] = sub[ct_col] - sub["Ct_ref"]
        results.append(sub)

    out = pd.concat(results, ignore_index=True)
    keep = [c for c in [sample_col, group_col, gene_col, ct_col, "Ct_ref", "dCt"] if c in out.columns]
    return out[keep]

def compute_ddct(df, reference_gene, control_group, target_genes=None,
                  sample_col="Sample", group_col="Group", gene_col="Gene", ct_col="Ct"):
    """
    Compute delta-Ct, delta-delta-Ct, and fold change (2^-ddCt) per sample.

    Parameters
    ----------
    df : DataFrame in long format (Sample, Group, Gene, Ct)
    reference_gene : str, name of housekeeping/reference gene
    control_group : str, name of the control group to normalize against
    target_genes : list of str or None. If None, uses all genes except reference_gene.

    Returns
    -------
    DataFrame with columns: Sample, Group, Gene, Ct, Ct_ref, dCt, ddCt, fold_change, log2_fold_change
    """
    df = df.copy()
    if target_genes is None:
        target_genes = [g for g in df[gene_col].unique() if g != reference_gene]

    # reference Ct per sample
    ref = df[df[gene_col] == reference_gene][[sample_col, ct_col]].rename(
        columns={ct_col: "Ct_ref"})
    if ref[sample_col].duplicated().any():
        raise ValueError(
            f"Multiple '{reference_gene}' Ct values found for the same sample. "
            "Average technical replicates before running this tool."
        )

    results = []
    for gene in target_genes:
        sub = df[df[gene_col] == gene].merge(ref, on=sample_col, how="left")
        if sub["Ct_ref"].isna().any():
            missing = sub[sub["Ct_ref"].isna()][sample_col].tolist()
            raise ValueError(f"No reference-gene Ct found for sample(s): {missing}")

        sub["dCt"] = sub[ct_col] - sub["Ct_ref"]

        # mean dCt of the control group for this gene defines the calibrator
        control_mean_dct = sub.loc[sub[group_col] == control_group, "dCt"].mean()
        if pd.isna(control_mean_dct):
            raise ValueError(f"Control group '{control_group}' not found for gene '{gene}'.")

        sub["ddCt"] = sub["dCt"] - control_mean_dct
        sub["fold_change"] = 2 ** (-sub["ddCt"])
        sub["log2_fold_change"] = np.log2(sub["fold_change"])
        results.append(sub)

    out = pd.concat(results, ignore_index=True)
    return out[[sample_col, group_col, gene_col, ct_col, "Ct_ref", "dCt", "ddCt",
                "fold_change", "log2_fold_change"]]


def run_stats(ddct_df, group_col="Group", gene_col="Gene", value_col="dCt",
              control_group=None, alpha=0.05, comparison_mode="pairwise"):
    """
    comparison_mode : "pairwise" (default) or "vs_control"
        "pairwise"   -- all group-vs-group comparisons (Tukey HSD, 3+ groups)
        "vs_control" -- each group compared only to control_group (Dunnett's
                        test, 3+ groups). Requires control_group to be set.
                        More statistically powerful than Tukey when you only
                        care about comparisons vs. control, since it isn't
                        spending correction on comparisons you're not making.
    """
    """
    Run stats per gene, on dCt values (NOT fold change -- fold change is
    log-normal, dCt is the appropriately-distributed value for parametric tests).

    Automatically chooses:
      - Welch's t-test if exactly 2 groups
      - One-way ANOVA + Tukey HSD post-hoc if 3+ groups

    Returns
    -------
    dict with keys:
      'summary'   : per-group per-gene mean fold change, SEM, n
      'pairwise'  : pairwise comparison p-values and significance stars
      'omnibus'   : ANOVA table per gene (only if 3+ groups)
    """
    genes = ddct_df[gene_col].unique()
    summary_rows = []
    pairwise_rows = []
    omnibus_rows = []

    for gene in genes:
        sub = ddct_df[ddct_df[gene_col] == gene]
        groups = sub[group_col].unique().tolist()
        n_groups = len(groups)

        # summary stats per group (on fold_change for reporting, dCt for testing)
        for g in groups:
            vals_fc = sub.loc[sub[group_col] == g, "fold_change"]
            summary_rows.append({
                "Gene": gene, "Group": g, "n": len(vals_fc),
                "mean_fold_change": vals_fc.mean(),
                "sem_fold_change": vals_fc.sem(),
                "sd_fold_change": vals_fc.std(),
            })

        if n_groups < 2:
            continue

        if n_groups == 2:
            g1, g2 = groups
            v1 = sub.loc[sub[group_col] == g1, value_col]
            v2 = sub.loc[sub[group_col] == g2, value_col]
            t_res = sp_stats.ttest_ind(v1, v2, equal_var=False)  # Welch's t-test
            p = t_res.pvalue
            pairwise_rows.append({
                "Gene": gene, "Group1": g1, "Group2": g2,
                "test": "Welch t-test", "p_value": p,
                "significance": p < alpha, "stars": _stars(p)
            })
        else:
            if comparison_mode == "vs_control" and control_group is None:
                raise ValueError("comparison_mode='vs_control' requires control_group to be set.")

            # one-way ANOVA (always computed, informational either way)
            model = ols(f"{value_col} ~ C({group_col})", data=sub).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            p_anova = anova_table["PR(>F)"].iloc[0]
            omnibus_rows.append({
                "Gene": gene, "test": "One-way ANOVA",
                "F": anova_table["F"].iloc[0], "p_value": p_anova,
                "significant": p_anova < alpha
            })

            if comparison_mode == "pairwise":
                # Tukey HSD post-hoc regardless of ANOVA significance (standard
                # practice to report, though usually only interpreted if omnibus p < alpha)
                tukey = pairwise_tukeyhsd(sub[value_col], sub[group_col], alpha=alpha)
                tukey_df = pd.DataFrame(data=tukey._results_table.data[1:],
                                         columns=tukey._results_table.data[0])
                for _, row in tukey_df.iterrows():
                    p = row["p-adj"]
                    pairwise_rows.append({
                        "Gene": gene, "Group1": row["group1"], "Group2": row["group2"],
                        "test": "Tukey HSD", "p_value": p,
                        "significant": bool(row["reject"]), "stars": _stars(p)
                    })

            elif comparison_mode == "vs_control":
                other_groups = [g for g in groups if g != control_group]
                control_vals = sub.loc[sub[group_col] == control_group, value_col].values
                other_samples = [sub.loc[sub[group_col] == g, value_col].values for g in other_groups]

                dunnett_res = dunnett(*other_samples, control=control_vals, alternative="two-sided")
                for g, p in zip(other_groups, dunnett_res.pvalue):
                    pairwise_rows.append({
                        "Gene": gene, "Group1": control_group, "Group2": g,
                        "test": "Dunnett", "p_value": p,
                        "significant": p < alpha, "stars": _stars(p)
                    })
            else:
                raise ValueError(f"Unknown comparison_mode: {comparison_mode}")

    return {
        "summary": pd.DataFrame(summary_rows),
        "pairwise": pd.DataFrame(pairwise_rows),
        "omnibus": pd.DataFrame(omnibus_rows) if omnibus_rows else None,
    }


def _stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


if __name__ == "__main__":
    # quick smoke test with synthetic data
    rng = np.random.default_rng(42)
    rows = []
    groups = {"Control": 18.0, "Treatment_Low": 17.5, "Treatment_High": 15.5}
    for group, target_mean in groups.items():
        for i in range(6):
            sample = f"{group}_{i+1}"
            ref_ct = rng.normal(18.2, 0.15)
            target_ct = rng.normal(target_mean, 0.2)
            rows.append({"Sample": sample, "Group": group, "Gene": "GAPDH", "Ct": ref_ct})
            rows.append({"Sample": sample, "Group": group, "Gene": "IL6", "Ct": target_ct})

    df = pd.DataFrame(rows)
    ddct = compute_ddct(df, reference_gene="GAPDH", control_group="Control")
    print(ddct.head(10))
    stats_out = run_stats(ddct, control_group="Control")
    print("\n--- Summary ---")
    print(stats_out["summary"])
    print("\n--- Omnibus (ANOVA) ---")
    print(stats_out["omnibus"])
    print("\n--- Pairwise (Tukey) ---")
    print(stats_out["pairwise"])
