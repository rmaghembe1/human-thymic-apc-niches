#!/usr/bin/env python3

from pathlib import Path
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
HVG_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_umap.h5ad"
PER_SAMPLE_DIR = PROJECT_DIR / "data/processed/heimli_filtered_per_sample"

OUT_RNA_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_RNA_only_for_markers.h5ad"
TAB_DIR = PROJECT_DIR / "results/tables"
TAB_DIR.mkdir(parents=True, exist_ok=True)

print("Reading downsampled UMAP object...")
hvg = sc.read_h5ad(HVG_H5AD)
selected_cells = set(hvg.obs_names)

pieces = []

print("Reconstructing selected cells from per-sample raw-count h5ad files...")
for f in sorted(PER_SAMPLE_DIR.glob("*_filtered.h5ad")):
    a = sc.read_h5ad(f)
    keep = [x for x in a.obs_names if x in selected_cells]

    if len(keep) == 0:
        continue

    a = a[keep, :].copy()

    for col in ["sample_prefix", "donor", "fraction", "data_type"]:
        a.obs[col] = a.obs[col].astype(str)

    print(f" - {f.name}: kept {a.n_obs} cells")
    pieces.append(a)
    gc.collect()

print("Concatenating selected cells...")
adata = ad.concat(
    pieces,
    join="inner",
    merge="same",
    index_unique=None
)

del pieces
gc.collect()

adata = adata[hvg.obs_names, :].copy()
adata.obsm["X_umap"] = hvg.obsm["X_umap"]
adata.obs["leiden_0_8"] = hvg.obs["leiden_0_8"].astype(str)

if "feature_type" not in adata.var.columns:
    raise ValueError("feature_type column is missing from adata.var. Cannot separate RNA and ADT features.")

print()
print("Feature-type counts:")
print(adata.var["feature_type"].value_counts().to_string())

rna_mask = adata.var["feature_type"].astype(str).eq("Gene Expression")
adt_mask = adata.var["feature_type"].astype(str).eq("Antibody Capture")

print()
print(f"RNA features: {int(rna_mask.sum())}")
print(f"ADT/protein features: {int(adt_mask.sum())}")

if rna_mask.sum() == 0:
    raise ValueError("No Gene Expression features found.")
if adt_mask.sum() == 0:
    print("Warning: no Antibody Capture features found. ADT table will be empty.")

print()
print("Building RNA-only object...")
rna = adata[:, rna_mask].copy()

print("Normalising and log-transforming RNA only...")
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

print("Saving RNA-only object...")
rna.write_h5ad(OUT_RNA_H5AD)

print("Running RNA-only rank_genes_groups...")
sc.tl.rank_genes_groups(
    rna,
    groupby="leiden_0_8",
    method="wilcoxon",
    n_genes=80
)

rna_rows = []
groups = rna.uns["rank_genes_groups"]["names"].dtype.names

for group in groups:
    df = sc.get.rank_genes_groups_df(rna, group=group).head(80).copy()
    df.insert(0, "cluster", group)
    df.insert(1, "rank", range(1, len(df) + 1))
    rna_rows.append(df)

rna_out = pd.concat(rna_rows, axis=0, ignore_index=True)
rna_file = TAB_DIR / "heimli_downsampled_rank_genes_by_leiden_RNA_only.tsv"
rna_out.to_csv(rna_file, sep="\t", index=False)

print("Building ADT/protein cluster mean table...")
if adt_mask.sum() > 0:
    adt = adata[:, adt_mask].copy()

    # Simple log1p transform for annotation-level ADT inspection.
    adt.X = adt.X.copy()
    adt.X.data = np.log1p(adt.X.data)

    adt_df = pd.DataFrame(
        adt.X.toarray() if hasattr(adt.X, "toarray") else adt.X,
        index=adt.obs_names,
        columns=adt.var_names
    )
    adt_df["leiden_0_8"] = adt.obs["leiden_0_8"].values

    adt_mean = adt_df.groupby("leiden_0_8", observed=True).mean().reset_index()

    adt_long = (
        adt_mean
        .melt(id_vars="leiden_0_8", var_name="ADT_feature", value_name="mean_log1p_ADT")
        .sort_values(["leiden_0_8", "mean_log1p_ADT"], ascending=[True, False])
    )

    adt_file = TAB_DIR / "heimli_downsampled_ADT_mean_by_leiden.tsv"
    adt_long.to_csv(adt_file, sep="\t", index=False)

    top_adt = adt_long.groupby("leiden_0_8").head(10)
else:
    adt_file = TAB_DIR / "heimli_downsampled_ADT_mean_by_leiden.tsv"
    pd.DataFrame(columns=["leiden_0_8", "ADT_feature", "mean_log1p_ADT"]).to_csv(adt_file, sep="\t", index=False)
    top_adt = pd.DataFrame()

print("Done.")
print(f"Wrote RNA-only object: {OUT_RNA_H5AD}")
print(f"Wrote RNA-only ranked genes: {rna_file}")
print(f"Wrote ADT/protein cluster means: {adt_file}")
print()
print("RNA-only top markers:")
print(rna_out.groupby("cluster").head(10).to_string(index=False))
print()
print("Top ADT/protein features:")
if len(top_adt) > 0:
    print(top_adt.to_string(index=False))
else:
    print("No ADT/protein features found.")
