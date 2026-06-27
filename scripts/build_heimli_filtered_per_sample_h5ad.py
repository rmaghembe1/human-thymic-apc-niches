#!/usr/bin/env python3

from pathlib import Path
import gzip
import re
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.io import mmread

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DATA_DIR = PROJECT_DIR / "data/raw/gse207206_heimli/extracted"
OUT_DIR = PROJECT_DIR / "data/processed/heimli_filtered_per_sample"
META_DIR = PROJECT_DIR / "metadata"
OUT_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

MIN_UMI = 500
MIN_GENES = 200
MAX_PCT_MT = 25

KEEP_FRACTIONS = {"APCenriched", "CD45depleted", "unenriched"}

def read_tsv_gz(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle]

def parse_prefix(prefix):
    m_donor = re.search(r"_(D[0-9]+)_", prefix)
    donor = m_donor.group(1) if m_donor else ""

    if "APCenriched" in prefix:
        fraction = "APCenriched"
    elif "CD45depleted" in prefix:
        fraction = "CD45depleted"
    elif "unenriched" in prefix:
        fraction = "unenriched"
    elif "PBMC" in prefix:
        fraction = "PBMC"
    else:
        fraction = "spatial_or_other"

    return donor, fraction

def load_features(feature_file):
    rows = read_tsv_gz(feature_file)
    gene_ids = []
    gene_symbols = []
    feature_types = []

    for row in rows:
        if len(row) >= 3:
            gene_ids.append(row[0])
            gene_symbols.append(row[1])
            feature_types.append(row[2])
        elif len(row) >= 2:
            gene_ids.append(row[0])
            gene_symbols.append(row[1])
            feature_types.append("Gene Expression")
        else:
            gene_ids.append(row[0])
            gene_symbols.append(row[0])
            feature_types.append("Gene Expression")

    return gene_ids, gene_symbols, feature_types

matrix_files = sorted(DATA_DIR.glob("*_raw_feature_bc_matrix_matrix.mtx.gz"))

summary_rows = []

for matrix_file in matrix_files:
    prefix = matrix_file.name.replace("_matrix.mtx.gz", "")
    donor, fraction = parse_prefix(prefix)

    if fraction not in KEEP_FRACTIONS:
        continue

    print(f"Processing {prefix} | donor={donor} | fraction={fraction}")

    barcode_file = DATA_DIR / f"{prefix}_barcodes.tsv.gz"
    feature_file = DATA_DIR / f"{prefix}_features.tsv.gz"

    if not barcode_file.exists() or not feature_file.exists():
        print(f"  Skipping: missing barcode or feature file")
        continue

    barcodes = [x[0] for x in read_tsv_gz(barcode_file)]
    gene_ids, gene_symbols, feature_types = load_features(feature_file)

    print("  Reading matrix...")
    mat = mmread(str(matrix_file)).tocsc()  # genes x barcodes, efficient for barcode filtering

    total_counts = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()

    keep_pre = (total_counts >= MIN_UMI) & (n_genes >= MIN_GENES)
    n_pre = int(keep_pre.sum())

    print(f"  Pre-filter kept {n_pre} of {mat.shape[1]} raw barcodes")

    if n_pre == 0:
        summary_rows.append({
            "prefix": prefix,
            "donor": donor,
            "fraction": fraction,
            "raw_barcodes": mat.shape[1],
            "pre_filter_barcodes": 0,
            "final_filtered_cells": 0,
            "median_umi_final": np.nan,
            "median_genes_final": np.nan,
            "median_pct_mt_final": np.nan,
            "output_h5ad": "",
        })
        del mat
        gc.collect()
        continue

    kept_barcodes = [barcodes[i] for i in np.where(keep_pre)[0]]
    X = mat[:, keep_pre].T.tocsr()  # cells x genes

    obs = pd.DataFrame(index=[f"{prefix}_{bc}" for bc in kept_barcodes])
    obs["dataset_id"] = "heimli_gse207206"
    obs["sample_prefix"] = prefix
    obs["donor"] = donor
    obs["fraction"] = fraction
    obs["data_type"] = "single_cell"
    obs["pre_filter_total_counts"] = total_counts[keep_pre]
    obs["pre_filter_n_genes"] = n_genes[keep_pre]

    var = pd.DataFrame(index=gene_symbols)
    var["gene_id"] = gene_ids
    var["gene_symbol"] = gene_symbols
    var["feature_type"] = feature_types

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

    print("  Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT].copy()

    out_file = OUT_DIR / f"{prefix}_filtered.h5ad"
    adata.write_h5ad(out_file)

    summary_rows.append({
        "prefix": prefix,
        "donor": donor,
        "fraction": fraction,
        "raw_barcodes": mat.shape[1],
        "pre_filter_barcodes": n_pre,
        "final_filtered_cells": adata.n_obs,
        "median_umi_final": round(float(adata.obs["total_counts"].median()), 2) if adata.n_obs else np.nan,
        "median_genes_final": round(float(adata.obs["n_genes_by_counts"].median()), 2) if adata.n_obs else np.nan,
        "median_pct_mt_final": round(float(adata.obs["pct_counts_mt"].median()), 3) if adata.n_obs else np.nan,
        "output_h5ad": str(out_file),
    })

    print(f"  Wrote {out_file} with {adata.n_obs} cells")

    del mat, X, adata
    gc.collect()

summary = pd.DataFrame(summary_rows)
summary_out = META_DIR / "heimli_gse207206_filtered_per_sample_qc_summary.tsv"
summary.to_csv(summary_out, sep="\t", index=False)

print("Done.")
print(f"Wrote summary: {summary_out}")
print()
print(summary.to_string(index=False))
