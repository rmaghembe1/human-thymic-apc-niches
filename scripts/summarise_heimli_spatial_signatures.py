#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR = PROJECT_DIR / "results/figures/heimli_spatial"
DOC_DIR = PROJECT_DIR / "docs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DOC_DIR.mkdir(parents=True, exist_ok=True)

sig_file = TAB_DIR / "heimli_spatial_signature_score_summary.tsv"
qc_file = PROJECT_DIR / "metadata/heimli_spatial_qc_summary.tsv"

sig = pd.read_csv(sig_file, sep="\t")
qc = pd.read_csv(qc_file, sep="\t")

mean_wide = sig.pivot(index="section_id", columns="signature", values="mean_score")
max_wide = sig.pivot(index="section_id", columns="signature", values="max_score")

mean_out = TAB_DIR / "heimli_spatial_signature_mean_score_matrix.tsv"
max_out = TAB_DIR / "heimli_spatial_signature_max_score_matrix.tsv"

mean_wide.to_csv(mean_out, sep="\t")
max_wide.to_csv(max_out, sep="\t")

ranked = (
    sig
    .sort_values(["section_id", "mean_score"], ascending=[True, False])
    .groupby("section_id", observed=True)
    .head(12)
    .copy()
)
ranked["rank_within_section"] = ranked.groupby("section_id", observed=True)["mean_score"].rank(
    method="first",
    ascending=False
).astype(int)

ranked_out = TAB_DIR / "heimli_spatial_signature_ranked_by_section.tsv"
ranked.to_csv(ranked_out, sep="\t", index=False)

# Simple derived axis scores.
axis_defs = {
    "axis_epithelial_TEC": ["sig_TEC", "sig_cTEC_like", "sig_mTEC_like"],
    "axis_APC_immune": ["sig_MHCII_AP", "sig_B_cell_APC", "sig_macrophage_APC", "sig_activated_DC", "sig_pDC_APC"],
    "axis_structural_stromal": ["sig_fibroblast", "sig_endothelial", "sig_perivascular"],
    "axis_T_cell_thymocyte": ["sig_T_cell_thymocyte"],
}

axis_rows = []
for section_id, sub in sig.groupby("section_id", observed=True):
    d = dict(zip(sub["signature"], sub["mean_score"]))
    row = {"section_id": section_id}
    for axis, members in axis_defs.items():
        vals = [d[m] for m in members if m in d]
        row[axis] = sum(vals) / len(vals) if vals else None
    axis_rows.append(row)

axis = pd.DataFrame(axis_rows)
axis_out = TAB_DIR / "heimli_spatial_niche_axis_scores_by_section.tsv"
axis.to_csv(axis_out, sep="\t", index=False)

def plot_heatmap(df, title, out_file):
    plot_df = df.copy()
    plot_df = plot_df.loc[sorted(plot_df.index), sorted(plot_df.columns)]

    fig_width = max(8, 0.55 * len(plot_df.columns))
    fig_height = max(4, 0.45 * len(plot_df.index))

    plt.figure(figsize=(fig_width, fig_height))
    plt.imshow(plot_df.values, aspect="auto")
    plt.xticks(range(len(plot_df.columns)), plot_df.columns, rotation=90)
    plt.yticks(range(len(plot_df.index)), plot_df.index)
    plt.title(title)
    plt.colorbar(label="score")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

plot_heatmap(
    mean_wide,
    "Heimli spatial mean signature scores",
    FIG_DIR / "heimli_spatial_mean_signature_score_heatmap.png"
)

plot_heatmap(
    max_wide,
    "Heimli spatial max signature scores",
    FIG_DIR / "heimli_spatial_max_signature_score_heatmap.png"
)

axis_plot = axis.set_index("section_id")
plot_heatmap(
    axis_plot,
    "Heimli spatial niche-axis scores",
    FIG_DIR / "heimli_spatial_niche_axis_score_heatmap.png"
)

summary_md = DOC_DIR / "heimli_spatial_signature_summary.md"

lines = []
lines.append("# Heimli spatial signature summary\n")
lines.append("## Spatial QC summary\n")
lines.append(qc.to_markdown(index=False))
lines.append("\n## Mean signature score matrix\n")
lines.append(mean_wide.reset_index().to_markdown(index=False))
lines.append("\n## Top-ranked signatures by section\n")
lines.append(ranked.sort_values(["section_id", "rank_within_section"]).to_markdown(index=False))
lines.append("\n## Niche-axis scores\n")
lines.append(axis.to_markdown(index=False))
lines.append("\n## Interpretation boundary\n")
lines.append("- These are first-pass signature scores, not formal spatial deconvolution.")
lines.append("- Scores should be interpreted comparatively across sections and spatial regions.")
lines.append("- Histology-overlay and spot-neighbourhood analyses remain required before making final spatial niche claims.")
summary_md.write_text("\n".join(lines))

print("Done.")
print(f"Wrote mean signature matrix: {mean_out}")
print(f"Wrote max signature matrix: {max_out}")
print(f"Wrote ranked signature table: {ranked_out}")
print(f"Wrote niche-axis table: {axis_out}")
print(f"Wrote summary document: {summary_md}")
print()
print("Niche-axis scores by section:")
print(axis.to_string(index=False))
print()
print("Top ranked signatures by section:")
print(ranked.sort_values(["section_id", "rank_within_section"]).to_string(index=False))
