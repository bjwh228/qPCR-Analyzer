"""
qPCR Analyzer - Streamlit App (v2, with import + QC pipeline)
==============================================================
Upload raw Ct data (any common format) -> assign groups -> QC (missing
values + technical replicate agreement) -> ddCt -> stats -> customizable figure.

Run locally:  streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import io

import hashlib

from qpcr_import import load_qpcr_file
from qpcr_qc import handle_missing_values, collapse_replicates
from qpcr_core import compute_ddct, compute_dct, run_stats
from qpcr_plot import plot_gene_with_points, DEFAULT_STYLE

st.set_page_config(page_title="qPCR Analyzer", layout="wide")
st.title("qPCR Analyzer")
st.caption("Upload raw Ct data → QC → ddCt → stats → publication-ready figures.")

# ================================================================ 1. Upload & parse
with st.sidebar:
    st.header("1. Upload data")
    uploaded = st.file_uploader(
        "Ct data file (long/wide CSV, Excel, or common instrument export)",
        type=["csv", "txt", "xlsx", "xls"])
    use_demo = st.checkbox("Use demo dataset instead", value=uploaded is None)

if use_demo:
    rng = np.random.default_rng(42)
    rows = []
    groups = {"Control": 18.0, "Treatment_Low": 17.5, "Treatment_High": 15.5}
    for group, target_mean in groups.items():
        for i in range(6):
            sample = f"{group}_{i+1}"
            for rep in range(2):  # simulate technical duplicates
                rows.append({"Sample": sample, "Group": group, "Gene": "GAPDH",
                            "Ct": rng.normal(18.2, 0.1)})
                rows.append({"Sample": sample, "Group": group, "Gene": "IL6",
                            "Ct": rng.normal(target_mean, 0.15)})
    raw_df = pd.DataFrame(rows)
    detected_format = "demo (long, with technical duplicates)"
elif uploaded is not None:
    try:
        raw_df, detected_format = load_qpcr_file(uploaded, filename=uploaded.name)
    except ValueError as e:
        st.error(str(e))
        st.stop()
else:
    st.info("Upload a file or check 'Use demo dataset' in the sidebar to get started.")
    st.stop()

st.sidebar.success(f"Detected format: {detected_format}")

n_missing = raw_df["Ct"].isna().sum()
if n_missing:
    st.sidebar.warning(f"{n_missing} missing/undetermined Ct value(s) found in the file.")

# Detect whether the underlying data source changed (new upload, or toggled
# demo on/off) -- if so, wipe out any group-assignment state tied to the old data.
if use_demo:
    data_source_signature = "demo"
else:
    data_source_signature = hashlib.md5(uploaded.getvalue()).hexdigest()

if st.session_state.get("data_source_signature") != data_source_signature:
    for key in ["sample_groups", "show_group_editor", "warning_dismissed"]:
        st.session_state.pop(key, None)
    st.session_state.data_source_signature = data_source_signature

st.sidebar.markdown("---")
analysis_mode = st.sidebar.radio(
    "Analysis type",
    ["Full ddCt analysis (needs groups + control)", "Just QC + dCt (no groups needed)"],
)
skip_groups = analysis_mode.startswith("Just QC")
has_groups_already = (
    "Group" in raw_df.columns
    and not raw_df["Group"].isna().all()
    and (raw_df["Group"].astype(str).str.strip() != "").any()
)

if skip_groups and has_groups_already:
    st.info(
        "This file already has group labels in it. If you want fold-change, statistics, "
        "and figures, switch to **Full ddCt analysis** in the sidebar."
    )
elif not skip_groups and not has_groups_already:
    st.info(
        "No group labels were found in this file. If you just want QC'd dCt values without "
        "assigning groups, you can switch to **Just QC + dCt** in the sidebar instead."
    )
# ================================================================ 2. Assign groups (if needed)
if skip_groups:
    edited_df = raw_df
else:
    if "sample_groups" not in st.session_state:
        unique_samples = raw_df[["Sample"]].drop_duplicates().reset_index(drop=True)
        unique_samples["Select"] = False
        if "Group" in raw_df.columns:
            existing = raw_df[["Sample", "Group"]].drop_duplicates(subset="Sample")
            unique_samples = unique_samples.merge(existing, on="Sample", how="left")
            unique_samples["Group"] = unique_samples["Group"].fillna("")
        else:
            unique_samples["Group"] = ""
        st.session_state.sample_groups = unique_samples

    if "show_group_editor" not in st.session_state:
        st.session_state.show_group_editor = (
            st.session_state.sample_groups["Group"].astype(str).str.strip() == ""
        ).any()

    if st.session_state.show_group_editor:
        st.header("1. Review imported data & assign groups")

        if "warning_dismissed" not in st.session_state:
            st.session_state.warning_dismissed = False

        if not st.session_state.warning_dismissed:
            warn_col1, warn_col2 = st.columns([20, 1])
            with warn_col1:
                st.warning("Select samples below and assign a group to all of them at once.")
            with warn_col2:
                st.write("")
                if st.button("✕", key="dismiss_warning"):
                    st.session_state.warning_dismissed = True
                    st.rerun()

        with st.form("group_assignment_form"):
            select_all = st.checkbox("Select all samples")
            display_df = st.session_state.sample_groups.copy()
            if select_all:
                display_df["Select"] = True

            edited_samples = st.data_editor(
                display_df,
                column_config={
                    "Select": st.column_config.CheckboxColumn("Select", width="small"),
                    "Sample": st.column_config.TextColumn("Sample", disabled=True, width="medium"),
                    "Group": st.column_config.TextColumn("Group", width="medium"),
                },
                hide_index=True,
                use_container_width=True,
                key="group_editor",
            )

            bulk_value = st.text_input("Group name to apply to selected samples")
            apply_clicked = st.form_submit_button("Apply to selected")

        st.session_state.sample_groups = edited_samples

        if apply_clicked and bulk_value.strip():
            mask = st.session_state.sample_groups["Select"]
            st.session_state.sample_groups.loc[mask, "Group"] = bulk_value.strip()
            st.session_state.sample_groups["Select"] = False
            st.rerun()

    edited_df = raw_df.drop(columns=["Group"], errors="ignore").merge(
        st.session_state.sample_groups[["Sample", "Group"]], on="Sample", how="left"
    )

    if (edited_df["Group"].astype(str).str.strip() == "").any():
        st.error("Some samples still have no Group assigned. Fill in the Group column above to continue.")
        st.stop()

    if not st.session_state.show_group_editor:
        st.header("1. Sample groups")
        st.success("All samples have a group assigned.")
        if st.button("Edit groups"):
            st.session_state.show_group_editor = True
            st.rerun()
        st.dataframe(edited_df, use_container_width=True)
# ================================================================ 3. Missing value handling
st.sidebar.markdown("---")
st.sidebar.header("2. Missing / undetermined values")
missing_strategy_label = st.sidebar.selectbox(
    "How should missing Ct values be handled?",
    ["Exclude just that well (keep other replicates)",
     "Impute a fixed Ct value",
     "Exclude the entire Sample+Gene if any replicate is missing"],
)
missing_strategy_map = {
    "Exclude just that well (keep other replicates)": "exclude_well",
    "Impute a fixed Ct value": "impute_fixed_ct",
    "Exclude the entire Sample+Gene if any replicate is missing": "exclude_sample_gene",
}
missing_strategy = missing_strategy_map[missing_strategy_label]

impute_ct = 40.0
if missing_strategy == "impute_fixed_ct":
    impute_ct = st.sidebar.number_input(
        "Fixed Ct value to impute (typically 40, your assay's cutoff)",
        min_value=20.0, max_value=50.0, value=40.0, step=0.5)

# ================================================================ 4. Replicate QC
st.sidebar.markdown("---")
st.sidebar.header("3. Technical replicate QC")
cutoff_type_label = st.sidebar.radio(
    "Disagreement measured in:", ["Fold change", "Raw Ct difference"], horizontal=True)

if cutoff_type_label == "Fold change":
    fold_cutoff = st.sidebar.slider(
        "Flag replicates that differ by this fold-change or more", 1.0, 4.0, 2.0, 0.1,
        help="A cutoff of 2.0 flags replicate pairs whose Ct values differ by 1 cycle or more "
             "(1 Ct = 2-fold).")
    cutoff_ct_equivalent = np.log2(fold_cutoff)
else:
    cutoff_ct_equivalent = st.sidebar.slider(
        "Flag replicates that differ by this many Ct or more", 0.1, 3.0, 0.5, 0.1,
        help="A common core-facility threshold is 0.3-0.5 Ct.")
    fold_cutoff = 2 ** cutoff_ct_equivalent

disagreement_label = st.sidebar.selectbox(
    "What to do with flagged replicate sets?",
    ["Flag but keep (average all replicates, mark for review)",
     "Exclude that Sample+Gene entirely",
     "Drop the outlier replicate and re-average (needs 3+ replicates)"],
)
disagreement_map = {
    "Flag but keep (average all replicates, mark for review)": "flag",
    "Exclude that Sample+Gene entirely": "exclude",
    "Drop the outlier replicate and re-average (needs 3+ replicates)": "drop_outlier",
}
on_disagreement = disagreement_map[disagreement_label]

# ================================================================ Run QC pipeline
clean_df, missing_log = handle_missing_values(
    edited_df, ct_col="Ct", strategy=missing_strategy, impute_ct=impute_ct)

collapsed_df, qc_report = collapse_replicates(
    clean_df, ct_col="Ct", cutoff=cutoff_ct_equivalent, cutoff_type="ct",
    on_disagreement=on_disagreement)

n_flagged = qc_report["flagged"].sum() if len(qc_report) else 0
if n_flagged:
    st.warning(f"⚠️ {n_flagged} Sample+Gene replicate set(s) flagged for disagreement — see the QC Report tab.")

# ================================================================ 5. Analysis settings
if collapsed_df.empty:
    st.error("No data left after QC filtering — try relaxing your missing-value or replicate cutoff settings.")
    st.stop()

genes = sorted(collapsed_df["Gene"].unique().tolist())

st.sidebar.markdown("---")
st.sidebar.header("4. Analysis settings")
reference_gene = st.sidebar.selectbox("Reference / housekeeping gene", genes)

if skip_groups:
    # ---- QC + dCt only, no groups/stats/figures ----
    target_genes = [g for g in genes if g != reference_gene]
    dct_df = compute_dct(collapsed_df, reference_gene=reference_gene, target_genes=target_genes)

    tab_qc, tab_data = st.tabs(["🔍 QC Report", "📄 dCt data"])

    with tab_qc:
        st.subheader("Missing value handling log")
        if len(missing_log):
            st.dataframe(missing_log)
        else:
            st.caption("No missing/undetermined values found.")
        st.subheader("Technical replicate QC")
        if len(qc_report):
            flagged_only = qc_report[qc_report["flagged"]]
            if len(flagged_only):
                st.dataframe(flagged_only)
            else:
                st.success("No replicate sets exceeded the cutoff.")
            with st.expander("Full replicate QC log"):
                st.dataframe(qc_report)

    with tab_data:
        st.subheader("Collapsed data used (post-QC)")
        st.dataframe(collapsed_df)
        st.subheader("dCt table")
        st.dataframe(dct_df)
        csv_buf = io.StringIO()
        dct_df.to_csv(csv_buf, index=False)
        st.download_button("Download dCt table (.csv)", csv_buf.getvalue(), file_name="qpcr_dct.csv")

    st.stop()
groups_list = sorted(collapsed_df["Group"].dropna().unique().tolist())
target_genes = st.sidebar.multiselect("Target gene(s) to analyze",
                                       [g for g in genes if g != reference_gene],
                                       default=[g for g in genes if g != reference_gene])
control_group = st.sidebar.selectbox("Control group", groups_list)
group_order = st.sidebar.multiselect("Group display order", groups_list, default=groups_list)
alpha = st.sidebar.slider("Significance threshold (alpha)", 0.01, 0.10, 0.05, 0.01)
comparison_mode_label = st.sidebar.radio(
    "Group comparisons",
    ["All pairwise comparisons", "Only vs. control group"],
    horizontal=False,
)
comparison_mode = "pairwise" if comparison_mode_label == "All pairwise comparisons" else "vs_control"

st.sidebar.markdown("---")
st.sidebar.header("5. Figure style")
palette_choice = st.sidebar.selectbox("Color palette", ["Default", "Blues", "Colorblind-safe", "Custom"])
if palette_choice == "Blues":
    palette = ["#08519c", "#3182bd", "#6baed6", "#bdd7e7"]
elif palette_choice == "Colorblind-safe":
    palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
elif palette_choice == "Custom":
    custom_colors = st.sidebar.text_input("Comma-separated hex colors", "#4C72B0,#DD8452,#55A868")
    palette = [c.strip() for c in custom_colors.split(",")]
else:
    palette = DEFAULT_STYLE["palette"]

font_family = st.sidebar.text_input("Font family", "DejaVu Sans")
font_size = st.sidebar.slider("Font size", 8, 20, 12)
y_scale = st.sidebar.radio("Y-axis scale", ["linear", "log2"], horizontal=True)
show_points = st.sidebar.checkbox("Show individual data points", value=True)
show_sig = st.sidebar.checkbox("Show significance brackets", value=True)
fig_w = st.sidebar.slider("Figure width (in)", 3, 10, 5)
fig_h = st.sidebar.slider("Figure height (in)", 3, 10, 4)
output_format = st.sidebar.selectbox("Export format", ["png", "pdf", "svg"])

# ================================================================ Run ddCt + stats
try:
    ddct = compute_ddct(collapsed_df, reference_gene=reference_gene, control_group=control_group,
                         target_genes=target_genes)
    stats_out = run_stats(ddct, control_group=control_group, alpha=alpha, comparison_mode=comparison_mode)
except ValueError as e:
    st.error(str(e))
    st.stop()

style = {
    **DEFAULT_STYLE,
    "palette": palette, "font_family": font_family, "font_size": font_size,
    "y_scale": y_scale, "show_points": show_points, "show_sig_brackets": show_sig,
    "fig_width": fig_w, "fig_height": fig_h, "output_format": output_format,
}

# ================================================================ Display
tab1, tab2, tab3, tab4 = st.tabs(["Figures", "Statistics", "QC Report", "Raw / ddCt data"])

with tab1:
    for gene in target_genes:
        st.subheader(gene)
        fig = plot_gene_with_points(gene, ddct, stats_out["summary"], stats_out["pairwise"],
                                     group_order=group_order, style=style)
        st.pyplot(fig)

        buf = io.BytesIO()
        fig.savefig(buf, format=output_format, dpi=style["dpi"], bbox_inches="tight")
        st.download_button(f"Download {gene} figure (.{output_format})", buf.getvalue(),
                            file_name=f"{gene}_qpcr.{output_format}", key=f"dl_{gene}")

with tab2:
    st.subheader("Group summary")
    st.dataframe(stats_out["summary"])

    if stats_out["omnibus"] is not None:
        st.subheader("one-way ANOVA")
        st.dataframe(stats_out["omnibus"])

    st.subheader("Pairwise comparisons")
    st.dataframe(stats_out["pairwise"])

    csv_buf = io.StringIO()
    stats_out["pairwise"].to_csv(csv_buf, index=False)
    st.download_button("Download stats (.csv)", csv_buf.getvalue(), file_name="qpcr_stats.csv")

with tab3:
    st.subheader("Missing value handling log")
    if len(missing_log):
        st.dataframe(missing_log)
    else:
        st.caption("No missing/undetermined values were found.")

    st.subheader("Technical replicate QC")
    st.caption(
        f"Cutoff: {cutoff_type_label.lower()} of {fold_cutoff:.2f}-fold "
        f"({cutoff_ct_equivalent:.2f} Ct-equivalent) · Action on disagreement: {disagreement_label}"
    )
    if len(qc_report):
        flagged_only = qc_report[qc_report["flagged"]]
        if len(flagged_only):
            st.markdown("**Flagged replicate sets:**")
            st.dataframe(flagged_only)
        else:
            st.success("No replicate sets exceeded the disagreement cutoff.")
        with st.expander("Show full replicate QC log (all Sample+Gene pairs)"):
            st.dataframe(qc_report)

    csv_buf3 = io.StringIO()
    qc_report.to_csv(csv_buf3, index=False)
    st.download_button("Download QC report (.csv)", csv_buf3.getvalue(), file_name="qpcr_qc_report.csv")

with tab4:
    st.subheader("Collapsed data used for analysis (post-QC)")
    st.dataframe(collapsed_df)
    st.subheader("ddCt / fold-change table")
    st.dataframe(ddct)
    csv_buf2 = io.StringIO()
    ddct.to_csv(csv_buf2, index=False)
    st.download_button("Download ddCt table (.csv)", csv_buf2.getvalue(), file_name="qpcr_ddct.csv")

st.markdown("---")
st.caption(
    "Method: ddCt. Stats run on dCt values (Welch's t-test for 2 groups, "
    "one-way ANOVA + Tukey HSD for 3+ groups) since fold-change is log-distributed. "
    "Technical replicates are QC'd and collapsed before analysis; see the QC Report tab for details."
)
