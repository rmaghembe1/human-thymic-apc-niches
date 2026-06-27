#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
IN_DIR = PROJECT_DIR / "data/processed/heimli_spatial_sections"
OUT_DIR = PROJECT_DIR / "data/processed/heimli_spatial_sections"
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR = PROJECT_DIR / "results/figures/heimli_spatial/copatterns"
DOC_DIR = PROJECT_DIR / "docs"

TAB_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
DOC_DIR.mkdir(parents=True, exist_ok=True)

SIGS = [
    "sig_TEC",
    "sig_cTEC_like",
    "sig_mTEC_like",
    "sig_MHCII_AP",
    "sig_B_cell_APC",
    "sig_macrophage_APC",
    "sig_activated_DC",
    "sig_pDC_APC",
    "sig_fibroblast",
    "sig_endothelial",
    "sig_perivascular",
    "sig_T_cell_thymocyte",
]

AXES = {
    "axis_epithelial_TEC": ["sig_TEC", "sig_cTEC_like", "sig_mTEC_like"],
    "axis_APC_immune": ["sig_MHCII_AP", "sig_B_cell_APC", "sig_macrophage_APC", "sig_activated_DC", "sig_pDC_APC"],
    "axis_structural_stromal": ["sig_fibroblast", "sig_endothelial", "sig_perivascular"],
    "axis_T_cell_thymocyte": ["sig_T_cell_thymocyte"],
}

COHIGH_RULES = {
    "cohigh_cTEC_MHCII": ("sig_cTEC_like_z", "sig_MHCII_AP_z"),
    "cohigh_mTEC_MHCII": ("sig_mTEC_like_z", "sig_MHCII_AP_z"),
    "cohigh_BcellAPC_MHCII": ("sig_B_cell_APC_z", "sig_MHCII_AP_z"),
    "cohigh_macrophageAPC_MHCII": ("sig_macrophage_APC_z", "sig_MHCII_AP_z"),
    "cohigh_pDC_MHCII": ("sig_pDC_APC_z", "sig_MHCII_AP_z"),
    "cohigh_fibroblast_TEC": ("sig_fibroblast_z", "axis_epithelial_TEC_z"),
    "cohigh_endothelial_APC": ("sig_endothelial_z", "axis_APC_immune_z"),
    "cohigh_thymocyte_APC": ("sig_T_cell_thymocyte_z", "axis_APC_immune_z"),
}

Z_THRESHOLD = 0.5

def zscore(values):
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if sd == 0 or np.isnan(sd):
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd

def plot_spatial_continuous(adata, section_id, col, out_file):
    df = adata.obs
    x = df["pxl_col_in_fullres"].astype(float).values
    y = df["pxl_row_in_fullres"].astype(float).values
    c = df[col].astype(float).values

    plt.figure(figsize=(5, 5))
    sca = plt.scatter(x, y, c=c, s=8)
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{section_id}: {col}")
    plt.colorbar(sca, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

def plot_spatial_categorical(adata, section_id, col, out_file):
    df = adata.obs
    x = df["pxl_col_in_fullres"].astype(float).values
    y = df["pxl_row_in_fullres"].astype(float).values
    cats = pd.Categorical(df[col].astype(str))
    codes = cats.codes

    plt.figure(figsize=(5, 5))
    plt.scatter(x, y, c=codes, s=8)
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{section_id}: {col}")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

spot_rows = []
axis_rows = []
cohigh_rows = []
corr_rows = []

files = sorted(IN_DIR.glob("*_spatial_scored.h5ad"))
if not files:
    raise FileNotFoundError(f"No spatial scored h5ad files found in {IN_DIR}")

for f in files:
    adata = sc.read_h5ad(f)
    section_id = str(adata.obs["section_id"].iloc[0])
    print(f"Processing co-patterns for {section_id}")

    present_sigs = [s for s in SIGS if s in adata.obs.columns]

    for sig in present_sigs:
        adata.obs[f"{sig}_z"] = zscore(adata.obs[sig].values)

    for axis_name, members in AXES.items():
        z_members = [f"{m}_z" for m in members if f"{m}_z" in adata.obs.columns]
        if len(z_members) == 0:
            adata.obs[f"{axis_name}_z"] = 0.0
        else:
            adata.obs[f"{axis_name}_z"] = adata.obs[z_members].mean(axis=1)

    sig_z_cols = [f"{s}_z" for s in present_sigs]
    axis_z_cols = [f"{a}_z" for a in AXES.keys()]

    sig_z = adata.obs[sig_z_cols]
    axis_z = adata.obs[axis_z_cols]

    adata.obs["dominant_signature"] = sig_z.idxmax(axis=1).str.replace("_z", "", regex=False)
    adata.obs["dominant_signature_z"] = sig_z.max(axis=1)
    adata.obs.loc[adata.obs["dominant_signature_z"] < Z_THRESHOLD, "dominant_signature"] = "low_mixed"

    adata.obs["dominant_axis"] = axis_z.idxmax(axis=1).str.replace("_z", "", regex=False)
    adata.obs["dominant_axis_z"] = axis_z.max(axis=1)
    adata.obs.loc[adata.obs["dominant_axis_z"] < Z_THRESHOLD, "dominant_axis"] = "low_mixed"

    for rule_name, (a, b) in COHIGH_RULES.items():
        if a in adata.obs.columns and b in adata.obs.columns:
            adata.obs[rule_name] = (
                (adata.obs[a] >= Z_THRESHOLD)
                & (adata.obs[b] >= Z_THRESHOLD)
            ).astype(int)
        else:
            adata.obs[rule_name] = 0

    out_h5ad = OUT_DIR / f"{f.stem.replace('_spatial_scored', '')}_spatial_copatterns.h5ad"
    adata.write_h5ad(out_h5ad)

    # Composition summaries
    for axis, n in adata.obs["dominant_axis"].value_counts().items():
        axis_rows.append({
            "section_id": section_id,
            "dominant_axis": axis,
            "n_spots": int(n),
            "pct_spots": round(100 * int(n) / adata.n_obs, 2),
        })

    for sig, n in adata.obs["dominant_signature"].value_counts().items():
        spot_rows.append({
            "section_id": section_id,
            "dominant_signature": sig,
            "n_spots": int(n),
            "pct_spots": round(100 * int(n) / adata.n_obs, 2),
        })

    for rule_name in COHIGH_RULES.keys():
        n = int(adata.obs[rule_name].sum())
        cohigh_rows.append({
            "section_id": section_id,
            "cohigh_rule": rule_name,
            "n_spots": n,
            "pct_spots": round(100 * n / adata.n_obs, 2),
        })

    # Signature correlation within section
    corr = adata.obs[present_sigs].corr(method="spearman")
    for sig1 in corr.index:
        for sig2 in corr.columns:
            if sig1 < sig2:
                corr_rows.append({
                    "section_id": section_id,
                    "signature_1": sig1,
                    "signature_2": sig2,
                    "spearman_rho": round(float(corr.loc[sig1, sig2]), 4),
                })

    # Plots
    for col in [
        "axis_epithelial_TEC_z",
        "axis_APC_immune_z",
        "axis_structural_stromal_z",
        "axis_T_cell_thymocyte_z",
        "sig_MHCII_AP_z",
        "sig_B_cell_APC_z",
        "sig_macrophage_APC_z",
        "sig_cTEC_like_z",
        "sig_mTEC_like_z",
    ]:
        if col in adata.obs.columns:
            plot_spatial_continuous(
                adata,
                section_id,
                col,
                FIG_DIR / f"{section_id}_{col}.png"
            )

    plot_spatial_categorical(
        adata,
        section_id,
        "dominant_axis",
        FIG_DIR / f"{section_id}_dominant_axis.png"
    )

    plot_spatial_categorical(
        adata,
        section_id,
        "dominant_signature",
        FIG_DIR / f"{section_id}_dominant_signature.png"
    )

axis_summary = pd.DataFrame(axis_rows)
sig_summary = pd.DataFrame(spot_rows)
cohigh_summary = pd.DataFrame(cohigh_rows)
corr_summary = pd.DataFrame(corr_rows)

axis_out = TAB_DIR / "heimli_spatial_dominant_axis_by_section.tsv"
sig_out = TAB_DIR / "heimli_spatial_dominant_signature_by_section.tsv"
cohigh_out = TAB_DIR / "heimli_spatial_cohigh_rules_by_section.tsv"
corr_out = TAB_DIR / "heimli_spatial_signature_spearman_correlations.tsv"

axis_summary.to_csv(axis_out, sep="\t", index=False)
sig_summary.to_csv(sig_out, sep="\t", index=False)
cohigh_summary.to_csv(cohigh_out, sep="\t", index=False)
corr_summary.to_csv(corr_out, sep="\t", index=False)

# Most relevant correlations for manuscript-level interpretation.
focus_pairs = [
    ("sig_cTEC_like", "sig_MHCII_AP"),
    ("sig_mTEC_like", "sig_MHCII_AP"),
    ("sig_B_cell_APC", "sig_MHCII_AP"),
    ("sig_macrophage_APC", "sig_MHCII_AP"),
    ("sig_pDC_APC", "sig_MHCII_AP"),
    ("sig_T_cell_thymocyte", "sig_MHCII_AP"),
    ("sig_fibroblast", "sig_TEC"),
    ("sig_endothelial", "sig_MHCII_AP"),
]

focus_rows = []
for section_id, sub in corr_summary.groupby("section_id", observed=True):
    for a, b in focus_pairs:
        row = sub[
            ((sub["signature_1"] == a) & (sub["signature_2"] == b))
            | ((sub["signature_1"] == b) & (sub["signature_2"] == a))
        ]
        if len(row) > 0:
            focus_rows.append(row.iloc[0].to_dict())

focus_corr = pd.DataFrame(focus_rows)
focus_corr_out = TAB_DIR / "heimli_spatial_focus_signature_correlations.tsv"
focus_corr.to_csv(focus_corr_out, sep="\t", index=False)

print("Done.")
print(f"Wrote dominant axis table: {axis_out}")
print(f"Wrote dominant signature table: {sig_out}")
print(f"Wrote co-high table: {cohigh_out}")
print(f"Wrote all correlation table: {corr_out}")
print(f"Wrote focus correlation table: {focus_corr_out}")
print()
print("Dominant axis by section:")
print(axis_summary.sort_values(["section_id", "pct_spots"], ascending=[True, False]).to_string(index=False))
print()
print("Co-high rules by section:")
print(cohigh_summary.sort_values(["section_id", "pct_spots"], ascending=[True, False]).to_string(index=False))
print()
print("Focus signature correlations:")
print(focus_corr.sort_values(["section_id", "signature_1", "signature_2"]).to_string(index=False))
