#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
TAB_DIR = PROJECT_DIR / "results/tables"
SPATIAL_H5AD_DIR = PROJECT_DIR / "data/processed/heimli_spatial_sections"
OUT_DIR = PROJECT_DIR / "results/figures/manuscript_final_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Automatically find the S2_B1 copattern h5ad file.
matches = sorted(SPATIAL_H5AD_DIR.glob("*S2_B1*_spatial_copatterns.h5ad"))
if not matches:
    raise FileNotFoundError(
        f"Could not find an S2_B1 spatial copattern h5ad in {SPATIAL_H5AD_DIR}. "
        "Expected a file like GSM6281325_S2_B1_spatial_copatterns.h5ad"
    )
SPATIAL_H5AD = matches[0]
print(f"Using spatial file: {SPATIAL_H5AD}")

cohigh = pd.read_csv(TAB_DIR / "heimli_spatial_cohigh_summary_across_sections.tsv", sep="\t")
corr = pd.read_csv(TAB_DIR / "heimli_spatial_focus_signature_correlations.tsv", sep="\t")
adata = sc.read_h5ad(SPATIAL_H5AD)

label_map = {
    "cohigh_BcellAPC_MHCII": "B-cell APC + MHC-II",
    "cohigh_mTEC_MHCII": "mTEC-like + MHC-II",
    "cohigh_macrophageAPC_MHCII": "Macrophage APC + MHC-II",
    "cohigh_pDC_MHCII": "pDC-like + MHC-II",
    "cohigh_endothelial_APC": "Endothelial + APC",
    "cohigh_fibroblast_TEC": "Fibroblast + TEC",
    "cohigh_cTEC_MHCII": "cTEC-like + MHC-II",
    "cohigh_thymocyte_APC": "Thymocyte + APC",
}

pair_label_map = {
    "sig_B_cell_APC :: sig_MHCII_AP": "B-cell APC : MHC-II",
    "sig_MHCII_AP :: sig_mTEC_like": "mTEC-like : MHC-II",
    "sig_MHCII_AP :: sig_macrophage_APC": "Macrophage APC : MHC-II",
    "sig_MHCII_AP :: sig_pDC_APC": "pDC-like : MHC-II",
    "sig_MHCII_AP :: sig_endothelial": "Endothelial : MHC-II",
    "sig_TEC :: sig_fibroblast": "Fibroblast : TEC",
    "sig_MHCII_AP :: sig_cTEC_like": "cTEC-like : MHC-II",
    "sig_MHCII_AP :: sig_T_cell_thymocyte": "T-cell/thymocyte : MHC-II",
}

spatial_panels = [
    ("A", "MHC-II/APC", "sig_MHCII_AP_z"),
    ("B", "B-cell APC", "sig_B_cell_APC_z"),
    ("C", "mTEC-like", "sig_mTEC_like_z"),
    ("D", "cTEC-like", "sig_cTEC_like_z"),
    ("E", "T-cell/thymocyte", "axis_T_cell_thymocyte_z"),
]

required_cols = [x[2] for x in spatial_panels]
missing = [c for c in required_cols if c not in adata.obs.columns]
if missing:
    raise ValueError(f"Missing expected columns in {SPATIAL_H5AD}: {missing}")

# Prepare Panel F data.
cohigh_plot = cohigh.copy()
cohigh_plot["label"] = cohigh_plot["cohigh_rule"].map(label_map).fillna(cohigh_plot["cohigh_rule"])
cohigh_plot = cohigh_plot.sort_values("mean_pct_spots", ascending=True)

# Prepare Panel G data.
corr_plot = corr.copy()
corr_plot["pair"] = corr_plot["signature_1"] + " :: " + corr_plot["signature_2"]
corr_plot["pair_label"] = corr_plot["pair"].map(pair_label_map).fillna(corr_plot["pair"])

pair_order = [
    "B-cell APC : MHC-II",
    "mTEC-like : MHC-II",
    "Macrophage APC : MHC-II",
    "pDC-like : MHC-II",
    "Endothelial : MHC-II",
    "Fibroblast : TEC",
    "cTEC-like : MHC-II",
    "T-cell/thymocyte : MHC-II",
]

section_order = ["S1_A1", "S1_B1", "S1_C1", "S1_D1", "S2_A1", "S2_B1", "S2_C1", "S2_D1"]

heat = (
    corr_plot
    .pivot(index="pair_label", columns="section_id", values="spearman_rho")
    .reindex(pair_order)
    .reindex(columns=section_order)
)

# Spatial coordinates.
x = adata.obs["pxl_col_in_fullres"].astype(float).values
y = adata.obs["pxl_row_in_fullres"].astype(float).values

# Build final figure.
fig = plt.figure(figsize=(18, 8.2))
gs = GridSpec(
    nrows=2,
    ncols=10,
    figure=fig,
    height_ratios=[1.0, 0.95],
    wspace=0.7,
    hspace=0.5
)

# Shared z-score range for A-E.
vmin, vmax = -3, 3
last_scatter = None

for i, (letter, title, col) in enumerate(spatial_panels):
    ax = fig.add_subplot(gs[0, i*2:(i+1)*2])
    vals = adata.obs[col].astype(float).values
    last_scatter = ax.scatter(x, y, c=vals, s=7, vmin=vmin, vmax=vmax)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(f"{letter}. {title}", fontsize=12, loc="left")

# Shared z-score colorbar for A-E.
cax = fig.add_axes([0.92, 0.56, 0.012, 0.30])
cb = fig.colorbar(last_scatter, cax=cax)
cb.set_label("z-score", fontsize=9)
cb.ax.tick_params(labelsize=8)

# Panel F: co-high bar plot.
ax_f = fig.add_subplot(gs[1, 0:4])
ax_f.barh(cohigh_plot["label"], cohigh_plot["mean_pct_spots"])
ax_f.errorbar(
    cohigh_plot["mean_pct_spots"],
    cohigh_plot["label"],
    xerr=[
        cohigh_plot["mean_pct_spots"] - cohigh_plot["min_pct_spots"],
        cohigh_plot["max_pct_spots"] - cohigh_plot["mean_pct_spots"],
    ],
    fmt="none",
    capsize=3,
)
ax_f.set_xlabel("Mean co-high spots across sections (%; min-max)", fontsize=10)
ax_f.set_ylabel("")
ax_f.set_title("F. Recurrent spatial co-high patterns", fontsize=12, loc="left")
ax_f.tick_params(axis="both", labelsize=9)
ax_f.grid(axis="x", alpha=0.25)
ax_f.set_xlim(0, max(cohigh_plot["max_pct_spots"]) * 1.15)

# Panel G: correlation heatmap.
ax_g = fig.add_subplot(gs[1, 4:10])
im = ax_g.imshow(heat.values, aspect="auto", vmin=-1, vmax=1)

ax_g.set_xticks(np.arange(len(heat.columns)))
ax_g.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=9)
ax_g.set_yticks(np.arange(len(heat.index)))
ax_g.set_yticklabels(heat.index, fontsize=9)
ax_g.set_title("G. Spatial signature correlations", fontsize=12, loc="left")

for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        val = heat.iloc[i, j]
        if pd.notna(val):
            ax_g.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

cbar = fig.colorbar(im, ax=ax_g, fraction=0.035, pad=0.02)
cbar.set_label("Spearman rho", fontsize=9)
cbar.ax.tick_params(labelsize=8)

out_png = OUT_DIR / "figure3_final_spatial_copatterns.png"
out_svg = OUT_DIR / "figure3_final_spatial_copatterns.svg"
out_pdf = OUT_DIR / "figure3_final_spatial_copatterns.pdf"

plt.savefig(out_png, dpi=600, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print("Done.")
print(f"Wrote: {out_png}")
print(f"Wrote: {out_svg}")
print(f"Wrote: {out_pdf}")
