#!/usr/bin/env python3

from pathlib import Path
import scanpy as sc
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
IN_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_refined_scores.h5ad"
TAB_DIR = PROJECT_DIR / "results/tables"
TAB_DIR.mkdir(parents=True, exist_ok=True)

print("Reading refined scored object...")
adata = sc.read_h5ad(IN_H5AD)

print("Running rank_genes_groups by Leiden cluster...")
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_0_6",
    method="wilcoxon",
    pts=True,
    n_genes=50
)

rows = []
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names

for group in groups:
    names = result["names"][group]
    scores = result["scores"][group]
    logfc = result["logfoldchanges"][group]
    pvals_adj = result["pvals_adj"][group]
    pct_nz_group = result["pts"][group] if "pts" in result else [None] * len(names)

    for rank, gene in enumerate(names, start=1):
        rows.append({
            "cluster": group,
            "rank": rank,
            "gene": gene,
            "score": scores[rank - 1],
            "logfoldchange": logfc[rank - 1],
            "pvals_adj": pvals_adj[rank - 1],
            "pct_nonzero_in_cluster": pct_nz_group[rank - 1] if pct_nz_group is not None else None,
        })

out = pd.DataFrame(rows)
out_file = TAB_DIR / "bautista_downsampled_rank_genes_by_leiden.tsv"
out.to_csv(out_file, sep="\t", index=False)

print(f"Wrote: {out_file}")
print()
print(out.groupby("cluster").head(10).to_string(index=False))
