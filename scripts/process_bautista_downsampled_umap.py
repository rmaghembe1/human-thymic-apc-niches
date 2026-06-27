#!/usr/bin/env python3

from pathlib import Path
import gc
import numpy as np
import scanpy as sc

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
IN_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_filtered.h5ad"
OUT_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_umap.h5ad"
FIG_DIR = PROJECT_DIR / "results/figures"
TAB_DIR = PROJECT_DIR / "results/tables"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

RANDOM_SEED = 20260625
MAX_CELLS_PER_SAMPLE = 1500

marker_sets = {
    "score_broad_TEC": ["EPCAM", "KRT8", "KRT18"],
    "score_cTEC": ["FOXN1", "PSMB11", "PRSS16", "DLL4"],
    "score_mTEC": ["KRT5", "KRT14", "AIRE", "FEZF2", "CD80", "CD40", "CCL21"],
    "score_mimetic_TEC": ["POU2F3", "TRPM5", "NEUROD1", "MYOD1", "FOXA2"],
    "score_antigen_presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M", "CIITA"],
    "score_DC_APC": ["HLA-DRA", "CD74", "ITGAX", "CD1C", "CLEC9A", "XCR1"],
    "score_B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA", "CD40"],
    "score_macrophage_APC": ["LST1", "C1QA", "C1QB", "C1QC", "AIF1", "CD68", "FCGR3A"],
}

print("Reading filtered Bautista object...")
adata_full = sc.read_h5ad(IN_H5AD)

for col in ["sample_prefix", "sample_short", "age_group", "age_value"]:
    adata_full.obs[col] = adata_full.obs[col].astype(str)

print("Balanced downsampling...")
rng = np.random.default_rng(RANDOM_SEED)
selected = []

for sample in sorted(adata_full.obs["sample_prefix"].unique()):
    sample_cells = adata_full.obs_names[adata_full.obs["sample_prefix"] == sample].to_numpy()
    n_take = min(MAX_CELLS_PER_SAMPLE, len(sample_cells))
    chosen = rng.choice(sample_cells, size=n_take, replace=False)
    selected.extend(chosen)
    print(f" - {sample}: selected {n_take} of {len(sample_cells)}")

adata = adata_full[selected, :].copy()
del adata_full
gc.collect()

print(f"Downsampled object: {adata.n_obs} cells x {adata.n_vars} genes")

print("Normalising and log-transforming...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("Scoring marker programmes...")
available_gene_names = set(adata.var_names)

for score_name, genes in marker_sets.items():
    present = [g for g in genes if g in available_gene_names]
    adata.uns[f"{score_name}_genes_present"] = present
    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name)
    else:
        adata.obs[score_name] = 0.0

print("Selecting highly variable genes...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1500,
    flavor="seurat",
    batch_key="sample_prefix"
)

print("Subsetting to HVGs...")
adata = adata[:, adata.var["highly_variable"]].copy()
gc.collect()

print("Scaling without zero-centering...")
sc.pp.scale(adata, max_value=10, zero_center=False)

print("PCA, neighbours, UMAP, Leiden...")
sc.tl.pca(adata, n_comps=25, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=12, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_0_6")

print("Writing downsampled processed object...")
adata.write_h5ad(OUT_H5AD)

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

sample_summary = (
    adata.obs
    .groupby(["sample_short", "age_group", "age_value"], observed=True)
    .size()
    .rename("n_cells_downsampled")
    .reset_index()
)

cluster_out = TAB_DIR / "bautista_downsampled_leiden_marker_score_summary.tsv"
sample_out = TAB_DIR / "bautista_downsampled_sample_summary.tsv"

cluster_summary.to_csv(cluster_out, sep="\t", index=False)
sample_summary.to_csv(sample_out, sep="\t", index=False)

print("Plotting UMAPs...")
sc.pl.umap(
    adata,
    color=["leiden_0_6"],
    legend_loc="on data",
    frameon=False,
    save="_bautista_downsampled_leiden.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["age_group"],
    frameon=False,
    save="_bautista_downsampled_age_group.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["sample_short"],
    frameon=False,
    save="_bautista_downsampled_sample.png",
    show=False
)

sc.pl.umap(
    adata,
    color=[
        "score_broad_TEC",
        "score_cTEC",
        "score_mTEC",
        "score_mimetic_TEC",
        "score_antigen_presentation",
        "score_DC_APC",
        "score_B_cell",
        "score_macrophage_APC",
    ],
    frameon=False,
    save="_bautista_downsampled_marker_scores.png",
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
