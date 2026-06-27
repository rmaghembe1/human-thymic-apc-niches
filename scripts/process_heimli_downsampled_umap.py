#!/usr/bin/env python3

from pathlib import Path
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
IN_DIR = PROJECT_DIR / "data/processed/heimli_filtered_per_sample"
OUT_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_umap.h5ad"
FIG_DIR = PROJECT_DIR / "results/figures"
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

RANDOM_SEED = 20260625
MAX_CELLS_PER_LIBRARY = 1200

marker_sets = {
    "score_TEC": ["EPCAM", "KRT8", "KRT18"],
    "score_cTEC": ["FOXN1", "PSMB11", "PRSS16", "DLL4"],
    "score_mTEC": ["KRT5", "KRT14", "AIRE", "FEZF2", "CCL21"],
    "score_MHCII_AP": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"],
    "score_cDC": ["ITGAX", "CD1C", "CLEC10A", "FCER1A", "CLEC9A", "XCR1", "CCR7", "LAMP3"],
    "score_pDC": ["IL3RA", "LILRA4", "TCF4", "IRF7", "GZMB"],
    "score_B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "CD37"],
    "score_plasma": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHA1"],
    "score_macrophage": ["LST1", "C1QA", "C1QB", "C1QC", "AIF1", "CD68", "FCGR3A"],
    "score_T_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2"],
    "score_thymocyte": ["CD7", "DNTT", "RAG1", "RAG2", "CD4", "CD8A"],
    "score_fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "score_endothelial": ["PECAM1", "VWF", "KDR", "ENG", "PLVAP"],
    "score_cycling": ["MKI67", "TOP2A", "PCNA", "UBE2C"],
}

rng = np.random.default_rng(RANDOM_SEED)

h5ad_files = sorted(IN_DIR.glob("*_filtered.h5ad"))
if not h5ad_files:
    raise FileNotFoundError(f"No filtered h5ad files found in {IN_DIR}")

pieces = []
sample_rows = []

print("Reading and balanced-downsampling Heimli per-sample objects...")

for f in h5ad_files:
    a = sc.read_h5ad(f)

    for col in ["sample_prefix", "donor", "fraction", "data_type"]:
        a.obs[col] = a.obs[col].astype(str)

    n_total = a.n_obs
    n_take = min(MAX_CELLS_PER_LIBRARY, n_total)

    if n_take < n_total:
        chosen = rng.choice(a.obs_names.to_numpy(), size=n_take, replace=False)
        a = a[chosen, :].copy()
    else:
        a = a.copy()

    sample_prefix = a.obs["sample_prefix"].iloc[0]
    donor = a.obs["donor"].iloc[0]
    fraction = a.obs["fraction"].iloc[0]

    print(f" - {sample_prefix}: selected {a.n_obs} of {n_total}")

    sample_rows.append({
        "sample_prefix": sample_prefix,
        "donor": donor,
        "fraction": fraction,
        "n_total_filtered_cells": n_total,
        "n_downsampled_cells": a.n_obs,
    })

    pieces.append(a)
    gc.collect()

print("Concatenating downsampled libraries...")
adata = ad.concat(
    pieces,
    join="inner",
    merge="same",
    index_unique=None
)

del pieces
gc.collect()

print(f"Downsampled Heimli object: {adata.n_obs} cells x {adata.n_vars} features")

print("Normalising and log-transforming...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("Scoring marker programmes...")
available = set(adata.var_names)

for score_name, genes in marker_sets.items():
    present = [g for g in genes if g in available]
    adata.uns[f"{score_name}_genes_present"] = present
    print(f" - {score_name}: {len(present)}/{len(genes)} markers present")
    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name)
    else:
        adata.obs[score_name] = 0.0

print("Selecting HVGs...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat",
    batch_key="sample_prefix"
)

print("Subsetting to HVGs for memory-safe PCA/UMAP...")
adata = adata[:, adata.var["highly_variable"]].copy()
gc.collect()

print("Scaling without zero-centering...")
sc.pp.scale(adata, max_value=10, zero_center=False)

print("PCA, neighbours, UMAP, Leiden...")
sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=25)
sc.tl.umap(adata)
sc.tl.leiden(
    adata,
    resolution=0.8,
    key_added="leiden_0_8",
    flavor="igraph",
    directed=False,
    n_iterations=2
)

print("Writing downsampled Heimli object...")
adata.write_h5ad(OUT_H5AD)

score_cols = list(marker_sets.keys())

cluster_summary = (
    adata.obs
    .groupby("leiden_0_8", observed=True)
    .agg(
        n_cells=("leiden_0_8", "size"),
        dominant_fraction=("fraction", lambda x: x.value_counts().idxmax()),
        dominant_donor=("donor", lambda x: x.value_counts().idxmax()),
        **{f"{col}_mean": (col, "mean") for col in score_cols}
    )
    .reset_index()
)

sample_summary = pd.DataFrame(sample_rows)

cluster_out = TAB_DIR / "heimli_downsampled_leiden_marker_score_summary.tsv"
sample_out = TAB_DIR / "heimli_downsampled_sample_summary.tsv"

cluster_summary.to_csv(cluster_out, sep="\t", index=False)
sample_summary.to_csv(sample_out, sep="\t", index=False)

print("Plotting UMAPs...")
sc.pl.umap(
    adata,
    color=["leiden_0_8"],
    legend_loc="on data",
    frameon=False,
    save="_heimli_downsampled_leiden.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["fraction", "donor"],
    frameon=False,
    save="_heimli_downsampled_fraction_donor.png",
    show=False
)

sc.pl.umap(
    adata,
    color=[
        "score_MHCII_AP",
        "score_cDC",
        "score_pDC",
        "score_B_cell",
        "score_plasma",
        "score_macrophage",
        "score_T_cell",
        "score_thymocyte",
        "score_TEC",
        "score_fibroblast",
        "score_endothelial",
    ],
    frameon=False,
    save="_heimli_downsampled_marker_scores.png",
    show=False
)

print("Done.")
print(f"Wrote downsampled object: {OUT_H5AD}")
print(f"Wrote cluster summary: {cluster_out}")
print(f"Wrote sample summary: {sample_out}")
print()
print("Sample summary:")
print(sample_summary.to_string(index=False))
print()
print("Cluster marker-score summary:")
print(cluster_summary.to_string(index=False))
