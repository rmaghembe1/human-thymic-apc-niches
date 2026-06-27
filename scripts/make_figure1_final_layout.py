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

BAUTISTA_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_annotated.h5ad"
HEIMLI_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_annotated.h5ad"

bautista = sc.read_h5ad(BAUTISTA_H5AD)
heimli = sc.read_h5ad(HEIMLI_H5AD)

def choose_col(adata, options, required=True):
    for col in options:
        if col in adata.obs.columns:
            return col
    if required:
        raise ValueError(f"Could not find any of these columns: {options}")
    return None

broad_col_b = choose_col(bautista, ["broad_class"])
sample_col_b = choose_col(
    bautista,
    ["sample_short", "sample_id", "sample", "sample_label", "sample_name", "sample_prefix"]
)

broad_col_h = choose_col(heimli, ["broad_class"])
fraction_col_h = choose_col(heimli, ["fraction"])

# Biological sample order for Bautista.
BAUTISTA_SAMPLE_ORDER = [
    "F_19wks",
    "F_23wks_1",
    "F_23wks_2",
    "P_6d",
    "P_10m_1",
    "P_10m_2",
    "A_25y",
]

# Simplify Heimli broad classes for Figure 1.
def simplify_heimli_class(x):
    x = str(x)

    immune_apc = {
        "B_cell_APC",
        "Macrophage_APC",
        "Activated_DC",
        "pDC_APC",
        "Activated_APC_like",
    }

    tec_epi = {
        "TEC_cTEC_like",
        "TEC_mTEC_like",
        "Epithelial_MHC_mixed",
        "Epithelial_stromal_mixed",
    }

    stromal_vascular = {
        "Fibroblast",
        "Endothelial",
        "Perivascular",
    }

    tcell_thymocyte = {
        "T_cell",
        "Thymocyte",
        "Thymocyte_mixed",
        "Cycling_thymocyte",
        "Immature_thymocyte",
        "Lymphoid_mixed",
    }

    mixed = {
        "Mixed_low_confidence",
    }

    contaminant = {
        "Contaminant",
    }

    if x in immune_apc:
        return "Immune APC"
    if x in tec_epi:
        return "TEC/epithelial"
    if x in stromal_vascular:
        return "Stromal/vascular"
    if x in tcell_thymocyte:
        return "T cell/thymocyte"
    if x in mixed:
        return "Mixed/low-confidence"
    if x in contaminant:
        return "Contaminant"

    return "Other"

heimli.obs["figure1_simplified_class"] = heimli.obs[broad_col_h].map(simplify_heimli_class).astype("category")

def get_umap(adata):
    if "X_umap" not in adata.obsm:
        raise ValueError("X_umap not found in AnnData object.")
    return adata.obsm["X_umap"]

def category_colors(categories, cmap_name="tab20"):
    cats = list(pd.Categorical(categories).categories)
    cmap = plt.get_cmap(cmap_name)
    colors = {cat: cmap(i % cmap.N) for i, cat in enumerate(cats)}
    return colors

def plot_umap(ax, adata, col, title, legend=True, legend_fontsize=7, point_size=4, legend_title=None):
    xy = get_umap(adata)
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

    if legend:
        handles = [
            Line2D(
                [0], [0],
                marker="o",
                color="w",
                label=str(cat),
                markerfacecolor=colors[cat],
                markersize=5,
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

def plot_stacked_composition(
    ax,
    adata,
    group_col,
    class_col,
    title,
    group_order=None,
    max_classes=None,
    legend_title="Class"
):
    df = adata.obs[[group_col, class_col]].copy()
    df[group_col] = df[group_col].astype(str)
    df[class_col] = df[class_col].astype(str)

    counts = (
        df
        .groupby([group_col, class_col], observed=True)
        .size()
        .reset_index(name="n")
    )

    if max_classes is not None:
        class_totals = counts.groupby(class_col, observed=True)["n"].sum().sort_values(ascending=False)
        keep = list(class_totals.head(max_classes).index)
        counts[class_col] = np.where(counts[class_col].isin(keep), counts[class_col], "Other")

    counts = (
        counts
        .groupby([group_col, class_col], observed=True)["n"]
        .sum()
        .reset_index()
    )

    wide = counts.pivot(index=group_col, columns=class_col, values="n").fillna(0)

    if group_order is not None:
        present = [g for g in group_order if g in wide.index]
        remaining = [g for g in wide.index if g not in present]
        wide = wide.loc[present + sorted(remaining)]
    else:
        wide = wide.loc[sorted(wide.index)]

    wide_pct = wide.div(wide.sum(axis=1), axis=0) * 100

    colors = category_colors(pd.Categorical(wide_pct.columns), cmap_name="tab20")
    bottom = np.zeros(wide_pct.shape[0])
    x = np.arange(wide_pct.shape[0])

    for cls in wide_pct.columns:
        vals = wide_pct[cls].values
        ax.bar(x, vals, bottom=bottom, label=cls, color=colors[cls], width=0.78)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(wide_pct.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Cells (%)", fontsize=9)
    ax.set_title(title, loc="left", fontsize=12)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="y", alpha=0.25)

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=7,
        title=legend_title,
        title_fontsize=7,
    )

fig = plt.figure(figsize=(18, 9.5))
gs = GridSpec(
    nrows=2,
    ncols=3,
    figure=fig,
    wspace=0.82,
    hspace=0.58
)

# Row 1: Bautista
ax_a = fig.add_subplot(gs[0, 0])
plot_umap(
    ax_a,
    bautista,
    broad_col_b,
    "A. Bautista stromal/TEC broad classes",
    legend=True,
    legend_fontsize=7,
    point_size=4
)

ax_b = fig.add_subplot(gs[0, 1])
plot_stacked_composition(
    ax_b,
    bautista,
    sample_col_b,
    broad_col_b,
    "B. Bautista sample composition",
    group_order=BAUTISTA_SAMPLE_ORDER,
    max_classes=8,
    legend_title="Class"
)

ax_c = fig.add_subplot(gs[0, 2])
plot_umap(
    ax_c,
    bautista,
    sample_col_b,
    "C. Bautista sample projection",
    legend=True,
    legend_fontsize=7,
    point_size=4,
    legend_title="Sample"
)

# Row 2: Heimli
ax_d = fig.add_subplot(gs[1, 0])
plot_umap(
    ax_d,
    heimli,
    "figure1_simplified_class",
    "D. Heimli simplified niche classes",
    legend=True,
    legend_fontsize=8,
    point_size=3,
    legend_title="Class"
)

ax_e = fig.add_subplot(gs[1, 1])
plot_stacked_composition(
    ax_e,
    heimli,
    fraction_col_h,
    "figure1_simplified_class",
    "E. Heimli fraction composition",
    group_order=["APCenriched", "CD45depleted", "unenriched"],
    max_classes=None,
    legend_title="Class"
)

ax_f = fig.add_subplot(gs[1, 2])
plot_umap(
    ax_f,
    heimli,
    fraction_col_h,
    "F. Heimli fraction projection",
    legend=True,
    legend_fontsize=8,
    point_size=3,
    legend_title="Fraction"
)

out_png = OUT_DIR / "figure1_final_dataset_architecture.png"
out_svg = OUT_DIR / "figure1_final_dataset_architecture.svg"
out_pdf = OUT_DIR / "figure1_final_dataset_architecture.pdf"

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
print(f"Bautista broad class: {broad_col_b}")
print(f"Bautista sample: {sample_col_b}")
print(f"Heimli broad class: {broad_col_h}")
print(f"Heimli fraction: {fraction_col_h}")
print()
print("Heimli simplified class counts:")
print(heimli.obs['figure1_simplified_class'].value_counts().to_string())
