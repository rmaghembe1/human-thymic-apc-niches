#!/usr/bin/env python3

from pathlib import Path
import scanpy as sc
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
FULL_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_filtered.h5ad"
UMAP_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_umap.h5ad"
OUT_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_refined_scores.h5ad"
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR = PROJECT_DIR / "results/figures"
TAB_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

marker_sets = {
    "refined_broad_TEC": ["EPCAM", "KRT8", "KRT18"],
    "refined_cTEC": ["FOXN1", "PSMB11", "PRSS16", "DLL4"],
    "refined_mTEC": ["KRT5", "KRT14", "AIRE", "FEZF2", "CCL21"],
    "refined_mimetic_TEC": ["POU2F3", "TRPM5", "NEUROD1", "MYOD1", "FOXA2"],
    "refined_MHCII_AP": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"],
    "refined_DC_specific": ["ITGAX", "CD1C", "CLEC9A", "XCR1", "FCER1A", "CCR7"],
    "refined_pDC_specific": ["IL3RA", "LILRA4", "TCF4", "IRF7", "GZMB"],
    "refined_B_specific": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "CD37"],
    "refined_macrophage_specific": ["LST1", "C1QA", "C1QB", "C1QC", "AIF1", "CD68", "FCGR3A"],
    "refined_T_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2"],
    "refined_fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "refined_endothelial": ["PECAM1", "VWF", "KDR", "ENG", "PLVAP"],
    "refined_cycling": ["MKI67", "TOP2A", "PCNA", "UBE2C"],
}

print("Reading downsampled UMAP object...")
umap = sc.read_h5ad(UMAP_H5AD)

print("Reading full filtered object and subsetting to UMAP cells...")
full = sc.read_h5ad(FULL_H5AD)
adata = full[umap.obs_names, :].copy()
del full

print("Normalising and log-transforming full-gene downsampled object...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("Copying UMAP coordinates and Leiden labels...")
adata.obsm["X_umap"] = umap.obsm["X_umap"]
adata.obs["leiden_0_6"] = umap.obs["leiden_0_6"].astype(str)

for col in ["sample_prefix", "sample_short", "age_group", "age_value"]:
    adata.obs[col] = adata.obs[col].astype(str)

print("Scoring refined marker programmes...")
available = set(adata.var_names)

for score_name, genes in marker_sets.items():
    present = [g for g in genes if g in available]
    adata.uns[f"{score_name}_genes_present"] = present
    print(f" - {score_name}: {len(present)}/{len(genes)} markers present")
    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name)
    else:
        adata.obs[score_name] = 0.0

score_cols = list(marker_sets.keys())

cluster_summary = (
    adata.obs
    .groupby("leiden_0_6", observed=True)
    .agg(
        n_cells=("leiden_0_6", "size"),
        dominant_age_group=("age_group", lambda x: x.value_counts().idxmax()),
        dominant_sample=("sample_short", lambda x: x.value_counts().idxmax()),
        **{f"{col}_mean": (col, "mean") for col in score_cols}
    )
    .reset_index()
)

out_table = TAB_DIR / "bautista_downsampled_refined_marker_score_summary.tsv"
cluster_summary.to_csv(out_table, sep="\t", index=False)

print("Plotting refined marker UMAPs...")
sc.pl.umap(
    adata,
    color=[
        "leiden_0_6",
        "refined_broad_TEC",
        "refined_cTEC",
        "refined_mTEC",
        "refined_mimetic_TEC",
        "refined_MHCII_AP",
        "refined_DC_specific",
        "refined_B_specific",
        "refined_macrophage_specific",
        "refined_T_cell",
        "refined_fibroblast",
        "refined_endothelial",
    ],
    frameon=False,
    save="_bautista_downsampled_refined_scores.png",
    show=False
)

adata.write_h5ad(OUT_H5AD)

print("Done.")
print(f"Wrote refined scored object: {OUT_H5AD}")
print(f"Wrote refined score summary: {out_table}")
print()
print(cluster_summary.to_string(index=False))
