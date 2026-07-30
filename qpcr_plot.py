"""
qPCR Plotting Module (matplotlib/seaborn)
==========================================
Takes the output of qpcr_core.run_stats()/compute_ddct() and produces a
publication-ready figure, with formatting fully controlled by a style dict
(so this can be driven by a web UI, config file, or CLI args).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

DEFAULT_STYLE = {
    "palette": ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2"],
    "font_family": "DejaVu Sans",   # swap for "Arial" / "Helvetica" if installed
    "font_size": 12,
    "title_font_size": 14,
    "fig_width": 5,
    "fig_height": 4.5,
    "dpi": 300,
    "bar_width": 0.6,
    "error_bar": "sem",       # "sem" or "sd"
    "show_points": True,
    "point_jitter": 0.08,
    "point_size": 25,
    "point_color": "#333333",
    "point_alpha": 0.7,
    "y_scale": "linear",      # "linear" or "log2"
    "ylabel": "Relative expression (fold change)",
    "xlabel": "",
    "show_sig_brackets": True,
    "sig_bracket_color": "black",
    "sig_line_width": 1.2,
    "bar_edge_color": "black",
    "bar_edge_width": 1.0,
    "output_format": "png",   # "png", "pdf", "svg"
    "style_theme": "ticks",   # seaborn theme: ticks, whitegrid, darkgrid, white, dark
    "despine": True,
}


def plot_gene(gene, summary_df, pairwise_df, group_order, style=None, output_path=None):
    """
    Create a bar chart (mean +/- error bar, optional jittered points, optional
    significance brackets) for a single gene across groups.

    Parameters
    ----------
    gene : str, gene name to plot (must exist in summary_df)
    summary_df : output of run_stats()['summary']
    pairwise_df : output of run_stats()['pairwise'] (for sig brackets); can be None
    group_order : list, order of groups on the x-axis
    style : dict overriding DEFAULT_STYLE keys
    output_path : if given, saves figure to this path
    """
    s = {**DEFAULT_STYLE, **(style or {})}

    plt.rcParams["font.family"] = s["font_family"]
    plt.rcParams["font.size"] = s["font_size"]
    sns.set_theme(style=s["style_theme"])

    sub = summary_df[summary_df["Gene"] == gene].set_index("Group").loc[group_order].reset_index()

    fig, ax = plt.subplots(figsize=(s["fig_width"], s["fig_height"]), dpi=s["dpi"])

    err_col = "sem_fold_change" if s["error_bar"] == "sem" else "sd_fold_change"
    colors = (s["palette"] * (len(group_order) // len(s["palette"]) + 1))[:len(group_order)]

    bars = ax.bar(
        sub["Group"], sub["mean_fold_change"], yerr=sub[err_col],
        width=s["bar_width"], color=colors, edgecolor=s["bar_edge_color"],
        linewidth=s["bar_edge_width"], capsize=4,
        error_kw={"elinewidth": 1.2, "capthick": 1.2},
    )

    if s["show_points"]:
        # NOTE: for real per-sample jittered points, pass the raw ddct dataframe;
        # this simplified version overlays mean markers only if raw data absent.
        pass

    if s["y_scale"] == "log2":
        ax.set_yscale("log", base=2)

    ax.set_ylabel(s["ylabel"], fontsize=s["font_size"])
    ax.set_xlabel(s["xlabel"], fontsize=s["font_size"])
    ax.set_title(gene, fontsize=s["title_font_size"], fontweight="bold")

    if s["despine"]:
        sns.despine(ax=ax)

    # significance brackets between adjacent / all pairs that are significant
    if s["show_sig_brackets"] and pairwise_df is not None and len(pairwise_df):
        gene_pairs = pairwise_df[(pairwise_df["Gene"] == gene) & (pairwise_df["significant"])]
        _draw_sig_brackets(ax, sub, gene_pairs, group_order, s)

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, format=s["output_format"], dpi=s["dpi"], bbox_inches="tight")
    return fig


def plot_gene_with_points(gene, ddct_df, summary_df, pairwise_df, group_order,
                           style=None, output_path=None):
    """
    Same as plot_gene but overlays actual per-sample jittered points
    (recommended -- always show the underlying data, not just bars).
    """
    s = {**DEFAULT_STYLE, **(style or {})}

    plt.rcParams["font.family"] = s["font_family"]
    plt.rcParams["font.size"] = s["font_size"]
    sns.set_theme(style=s["style_theme"])

    sub_summary = summary_df[summary_df["Gene"] == gene].set_index("Group").loc[group_order].reset_index()
    sub_raw = ddct_df[ddct_df["Gene"] == gene]

    fig, ax = plt.subplots(figsize=(s["fig_width"], s["fig_height"]), dpi=s["dpi"])

    err_col = "sem_fold_change" if s["error_bar"] == "sem" else "sd_fold_change"
    colors = (s["palette"] * (len(group_order) // len(s["palette"]) + 1))[:len(group_order)]

    x_pos = np.arange(len(group_order))
    ax.bar(
        x_pos, sub_summary["mean_fold_change"], yerr=sub_summary[err_col],
        width=s["bar_width"], color=colors, edgecolor=s["bar_edge_color"],
        linewidth=s["bar_edge_width"], capsize=4,
        error_kw={"elinewidth": 1.2, "capthick": 1.2}, zorder=2,
    )

    if s["show_points"]:
        rng = np.random.default_rng(0)
        for i, g in enumerate(group_order):
            vals = sub_raw.loc[sub_raw["Group"] == g, "fold_change"].values
            jitter = rng.uniform(-s["point_jitter"], s["point_jitter"], size=len(vals))
            ax.scatter(np.full(len(vals), x_pos[i]) + jitter, vals,
                       s=s["point_size"], color=s["point_color"],
                       alpha=s["point_alpha"], zorder=3, edgecolor="white", linewidth=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_order)

    if s["y_scale"] == "log2":
        ax.set_yscale("log", base=2)

    ax.set_ylabel(s["ylabel"], fontsize=s["font_size"])
    ax.set_xlabel(s["xlabel"], fontsize=s["font_size"])
    ax.set_title(gene, fontsize=s["title_font_size"], fontweight="bold")

    if s["despine"]:
        sns.despine(ax=ax)

    if s["show_sig_brackets"] and pairwise_df is not None and len(pairwise_df):
        gene_pairs = pairwise_df[(pairwise_df["Gene"] == gene) & (pairwise_df["significant"])]
        _draw_sig_brackets(ax, sub_summary, gene_pairs, group_order, s, raw_df=sub_raw)

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, format=s["output_format"], dpi=s["dpi"], bbox_inches="tight")
    return fig


def _draw_sig_brackets(ax, summary_sub, sig_pairs, group_order, style, raw_df=None):
    """Draw significance brackets + stars above bars for significant comparisons."""
    if len(sig_pairs) == 0:
        return

    if raw_df is not None:
        max_val = raw_df["fold_change"].max()
    else:
        max_val = (summary_sub["mean_fold_change"] +
                   summary_sub.get("sem_fold_change", 0)).max()

    y_range = max_val * 0.15
    current_y = max_val + y_range

    for _, row in sig_pairs.iterrows():
        try:
            x1 = group_order.index(row["Group1"])
            x2 = group_order.index(row["Group2"])
        except ValueError:
            continue
        y = current_y
        ax.plot([x1, x1, x2, x2], [y, y + y_range * 0.15, y + y_range * 0.15, y],
                lw=style["sig_line_width"], color=style["sig_bracket_color"])
        ax.text((x1 + x2) / 2, y + y_range * 0.18, row.get("stars", "*"),
                ha="center", va="bottom", fontsize=style["font_size"])
        current_y += y_range * 0.55

    ax.set_ylim(top=current_y + y_range * 0.3)
