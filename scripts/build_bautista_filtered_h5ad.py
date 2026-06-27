#!/usr/bin/env python3

from pathlib import Path
import gzip
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DATA_DIR = PROJECT_DIR / "data/raw/gse147520_bautista/extracted"
OUT_DIR = PROJECT_DIR / "data/processed"
META_DIR = PROJECT_DIR / "metadata"
OUT_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

samples = [
    ("GSM4466780_F_19wks", "fetal", "19wks", "F_19wks"),
    ("GSM4466781_F_23wks_1", "fetal", "23wks", "F_23wks_1"),
    ("GSM4466782_F_23wks_2", "fetal", "23wks", "F_23wks_2"),
    ("GSM4466783_P_6d", "postnatal", "6d", "P_6d"),
    ("GSM4466784_P_10m_1", "postnatal", "10m", "P_10m_1"),
    ("GSM4466785_P_10m_2", "postnatal", "10m", "P_10m_2"),
    ("GSM4466786_A_25y", "adult", "25y", "A_25y"),
]

def read_tsv_gz(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle]

def find_feature_file(prefix):
    genes = DATA_DIR / f"{prefix}_genes.tsv.gz"
    features = DATA_DIR / f"{prefix}_features.tsv.gz"
    if genes.exists():
        return genes
    if features.exists():
        return features
    raise FileNotFoundError(f"No genes/features file found for {prefix}")

def load_sample(prefix, age_group, age_value, sample_short):
    matrix_file = DATA_DIR / f"{prefix}_matrix.mtx.gz"
    barcode_file = DATA_DIR / f"{prefix}_barcodes.tsv.gz"
    feature_file = find_feature_file(prefix)

    matrix = mmread(str(matrix_file)).tocsr().T  # cells x genes

    barcodes = [x[0] for x in read_tsv_gz(barcode_file)]

    feature_rows = read_tsv_gz(feature_file)
    gene_ids = []
    gene_symbols = []
    for row in feature_rows:
        if len(row) >= 2:
            gene_ids.append(row[0])
            gene_symbols.append(row[1])
        else:
            gene_ids.append(row[0])
            gene_symbols.append(row[0])

    obs = pd.DataFrame(index=[f"{prefix}_{bc}" for bc in barcodes])
    obs["sample_prefix"] = prefix
    obs["sample_short"] = sample_short
    obs["age_group"] = age_group
    obs["age_value"] = age_value
    obs["dataset_id"] = "bautista_gse147520"

    var = pd.DataFrame(index=gene_symbols)
    var["gene_id"] = gene_ids
    var["gene_symbol"] = gene_symbols

    # Make duplicated gene symbols unique while preserving symbol column
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata

print("Loading Bautista samples...")
adatas = []
for sample in samples:
    print(" -", sample[0])
    adatas.append(load_sample(*sample))

print("Concatenating...")
adata = ad.concat(adatas, join="outer", label="batch", keys=[s[0] for s in samples], index_unique=None)

# Store raw counts
adata.layers["counts"] = adata.X.copy()

print("Calculating QC metrics...")
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

raw_out = OUT_DIR / "bautista_gse147520_raw_merged.h5ad"
adata.write_h5ad(raw_out)

# First-pass feasibility QC, deliberately conservative
adata_f = adata.copy()
sc.pp.filter_cells(adata_f, min_genes=200)
sc.pp.filter_genes(adata_f, min_cells=3)

# Additional broad QC
adata_f = adata_f[adata_f.obs["pct_counts_mt"] < 20].copy()
adata_f = adata_f[adata_f.obs["total_counts"] >= 500].copy()

filtered_out = OUT_DIR / "bautista_gse147520_filtered.h5ad"
adata_f.write_h5ad(filtered_out)

qc_summary = (
    adata.obs.groupby(["sample_prefix", "age_group", "age_value"])
    .agg(
        raw_cells=("sample_prefix", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_umi=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .reset_index()
)

filtered_counts = (
    adata_f.obs.groupby(["sample_prefix"])
    .size()
    .rename("filtered_cells")
    .reset_index()
)

qc_summary = qc_summary.merge(filtered_counts, on="sample_prefix", how="left")
qc_summary["filtered_cells"] = qc_summary["filtered_cells"].fillna(0).astype(int)

qc_out = META_DIR / "bautista_gse147520_qc_summary.tsv"
qc_summary.to_csv(qc_out, sep="\t", index=False)

print("Done.")
print(f"Wrote raw merged object: {raw_out}")
print(f"Wrote filtered object: {filtered_out}")
print(f"Wrote QC summary: {qc_out}")
print()
print(qc_summary.to_string(index=False))
