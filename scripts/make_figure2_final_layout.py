#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
OUT_DIR = PROJECT_DIR / "results/figures/manuscript_final_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HEIMLI_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_annotated.h5ad"
adata = sc.read_h5ad(HEIMLI_H5AD)

def choose_col(adata, options, required=True):
    for col in options:
        if col in adata.obs.columns:
            return col
    if required:
        raise ValueError(f"Could not find any of these columns: {options}")
    return None

broad_col = choose_col(adata, ["broad_class"])
confidence_col = choose_col(adata, ["annotation_confidence"])

score_cols = {
    "MHC-II/APC": ["score_MHCII_AP", "refined_MHCII_AP", "sig_MHCII_AP"],
    "B-cell APC": ["score_B_cell", "refined_B_specific", "score_B_cell_APC", "sig_B_cell_APC"],
    "Macrophage/APC": ["score_macrophage", "refined_macrophage_specific", "score_macrophage_APC", "sig_macrophage_APC"],
    "pDC/APC": ["score_pDC", "refined_pDC_specific", "score_pDC_APC", "sig_pDC_APC"],
    "cTEC-like": ["score_cTEC", "refined_cTEC", "sig_cTEC_like"],
    "mTEC-like": ["score_mTEC", "refined_mTEC", "sig_mTEC_like"],
}

resolved_scores = {}
for label, options in score_cols.items():
    resolved_scores[label] = choose_col(adata, options)

if "X_umap" not in adata.obsm:
    raise ValueError("X_umap not found in AnnData object.")

xy = adata.obsm["X_umap"]

def category_colors(categories, cmap_name="tab20"):
    cats = list(pd.Categorical(categories).categories)
    cmap = plt.get_cmap(cmap_name)
    colors = {cat: cmap(i % cmap.N) for i, cat in enumerate(cats)}
    return colors

def plot_categorical_umap(ax, col, title, legend_fontsize=5.8, point_size=3, legend_title="Class"):
    vals = pd.Categorical(adata.obs[col].astype(str))
    colors = category_colors(vals)
    color_values = [colors[v] for v in vals]

    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=color_values,
        s=point_size,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(title, loc="left", fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")

    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=str(cat),
            markerfacecolor=colors[cat],
            markersize=4,
        )
        for cat in colors.keys()
    ]

    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=legend_fontsize,
        title=legend_title,
        title_fontsize=legend_fontsize,
        borderaxespad=0.0,
    )

def robust_limits(values):
    values = np.asarray(values, dtype=float)
    lo = np.nanpercentile(values, 1)
    hi = np.nanpercentile(values, 99)
    if lo == hi:
        lo, hi = np.nanmin(values), np.nanmax(values)
    return lo, hi

def plot_score_umap(ax, col, title, point_size=3):
    vals = adata.obs[col].astype(float).values
    vmin, vmax = robust_limits(vals)

    sca = ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=vals,
        s=point_size,
        linewidths=0,
        rasterized=True,
        vmin=vmin,
        vmax=vmax
    )

    ax.set_title(title, loc="left", fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")

    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

fig = plt.figure(figsize=(18, 9.5))
gs = GridSpec(
    nrows=2,
    ncols=4,
    figure=fig,
    wspace=0.78,
    hspace=0.45
)

# A: detailed annotation map
ax_a = fig.add_subplot(gs[0, 0])
plot_categorical_umap(
    ax_a,
    broad_col,
    "A. Heimli annotated broad classes",
    legend_fontsize=5.7,
    point_size=3,
    legend_title="Class"
)

# B-G marker scores
panel_defs = [
    ("B", "MHC-II/APC", resolved_scores["MHC-II/APC"]),
    ("C", "B-cell APC", resolved_scores["B-cell APC"]),
    ("D", "Macrophage/APC", resolved_scores["Macrophage/APC"]),
    ("E", "pDC/APC", resolved_scores["pDC/APC"]),
    ("F", "cTEC-like", resolved_scores["cTEC-like"]),
    ("G", "mTEC-like", resolved_scores["mTEC-like"]),
]

positions = [
    gs[0, 1],
    gs[0, 2],
    gs[0, 3],
    gs[1, 0],
    gs[1, 1],
    gs[1, 2],
]

for (letter, label, col), pos in zip(panel_defs, positions):
    ax = fig.add_subplot(pos)
    plot_score_umap(
        ax,
        col,
        f"{letter}. {label} score",
        point_size=3
    )

# H: annotation confidence
ax_h = fig.add_subplot(gs[1, 3])
plot_categorical_umap(
    ax_h,
    confidence_col,
    "H. Annotation confidence",
    legend_fontsize=8,
    point_size=3,
    legend_title="Confidence"
)

out_png = OUT_DIR / "figure2_final_marker_programmes.png"
out_svg = OUT_DIR / "figure2_final_marker_programmes.svg"
out_pdf = OUT_DIR / "figure2_final_marker_programmes.pdf"

plt.savefig(out_png, dpi=600, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print("Done.")
print(f"Wrote: {out_png}")
print(f"Wrote: {out_svg}")
print(f"Wrote: {out_pdf}")
print()
print("Columns used:")
print(f"Broad class: {broad_col}")
print(f"Annotation confidence: {confidence_col}")
for label, col in resolved_scores.items():
    print(f"{label}: {col}")
