#!/usr/bin/env python3

from pathlib import Path
import gzip
import re
import json
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.io import mmread
import matplotlib.pyplot as plt

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DATA_DIR = PROJECT_DIR / "data/raw/gse207206_heimli/extracted"
OUT_DIR = PROJECT_DIR / "data/processed/heimli_spatial_sections"
FIG_DIR = PROJECT_DIR / "results/figures/heimli_spatial"
TAB_DIR = PROJECT_DIR / "results/tables"
META_DIR = PROJECT_DIR / "metadata"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

MIN_TOTAL_COUNTS = 100
MIN_GENES = 50
MAX_PCT_MT = 30

signature_sets = {
    "sig_TEC": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "sig_cTEC_like": ["CCL25", "PRSS16", "PSMB11", "FOXN1", "CTSV"],
    "sig_mTEC_like": ["KRT5", "KRT14", "AIRE", "FEZF2", "CCL19", "CALML5", "S100A14"],
    "sig_MHCII_AP": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"],
    "sig_B_cell_APC": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"],
    "sig_macrophage_APC": ["LYZ", "CST3", "LST1", "C1QA", "C1QB", "C1QC", "AIF1", "CD68"],
    "sig_activated_DC": ["CCR7", "LAMP3", "FSCN1", "BIRC3", "TNFAIP2"],
    "sig_pDC_APC": ["LILRA4", "IL3RA", "GZMB", "PLD4", "IRF7", "JCHAIN"],
    "sig_fibroblast": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "PTGDS"],
    "sig_endothelial": ["PECAM1", "VWF", "RAMP2", "ENG", "TM4SF1", "SPARCL1"],
    "sig_perivascular": ["MYL9", "TAGLN", "ACTA2", "CALD1", "RGS5", "PDGFRB"],
    "sig_T_cell_thymocyte": ["CD3D", "CD3E", "TRAC", "CD7", "DNTT", "RAG1", "RAG2", "PTCRA"],
}

def read_tsv_gz(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle]

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

def load_positions(path):
    # Old Space Ranger format:
    # barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
    pos = pd.read_csv(path, header=None)
    pos.columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    pos["barcode"] = pos["barcode"].astype(str)
    return pos

def plot_spatial_score(adata, section_id, color_col, out_file):
    df = adata.obs.copy()
    x = df["pxl_col_in_fullres"].astype(float).values
    y = df["pxl_row_in_fullres"].astype(float).values
    c = df[color_col].astype(float).values

    plt.figure(figsize=(5, 5))
    sca = plt.scatter(x, y, c=c, s=8)
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"{section_id}: {color_col}")
    plt.colorbar(sca, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

matrix_files = sorted(DATA_DIR.glob("GSM62813*_S*_raw_feature_bc_matrix_matrix.mtx.gz"))

summary_rows = []
sig_rows = []

for matrix_file in matrix_files:
    prefix = matrix_file.name.replace("_raw_feature_bc_matrix_matrix.mtx.gz", "")

    # Spatial prefixes look like GSM6281320_S1_A1
    if not re.search(r"_S[12]_[ABCD]1$", prefix):
        continue

    section_id = prefix.split("_", 1)[1]

    barcode_file = DATA_DIR / f"{prefix}_raw_feature_bc_matrix_barcodes.tsv.gz"
    feature_file = DATA_DIR / f"{prefix}_raw_feature_bc_matrix_features.tsv.gz"
    positions_file = DATA_DIR / f"{prefix}_spatial_tissue_positions_list.csv.gz"
    scale_file = DATA_DIR / f"{prefix}_spatial_scalefactors_json.json.gz"

    print(f"Processing spatial section {section_id} ({prefix})")

    if not barcode_file.exists() or not feature_file.exists() or not positions_file.exists():
        print("  Skipping because required matrix/feature/barcode/position files are missing.")
        continue

    barcodes = [x[0] for x in read_tsv_gz(barcode_file)]
    gene_ids, gene_symbols, feature_types = load_features(feature_file)

    print("  Reading matrix...")
    mat = mmread(str(matrix_file)).T.tocsr()  # spots x genes

    obs = pd.DataFrame(index=barcodes)
    obs["barcode"] = barcodes
    obs["dataset_id"] = "heimli_gse207206"
    obs["data_type"] = "spatial_visium"
    obs["sample_prefix"] = prefix
    obs["section_id"] = section_id
    obs["slide"] = section_id.split("_")[0]
    obs["capture_area"] = section_id.split("_")[1]

    var = pd.DataFrame(index=gene_symbols)
    var["gene_id"] = gene_ids
    var["gene_symbol"] = gene_symbols
    var["feature_type"] = feature_types

    adata = ad.AnnData(X=mat, obs=obs, var=var)
    adata.var_names_make_unique()

    positions = load_positions(positions_file).set_index("barcode")
    adata.obs = adata.obs.join(positions, on="barcode")

    if adata.obs["in_tissue"].isna().any():
        n_missing = int(adata.obs["in_tissue"].isna().sum())
        print(f"  Warning: {n_missing} barcodes missing spatial position metadata")

    adata = adata[adata.obs["in_tissue"].fillna(0).astype(int) == 1, :].copy()

    # Keep RNA features only if feature_type exists and contains Gene Expression.
    if "feature_type" in adata.var.columns:
        ge_mask = adata.var["feature_type"].astype(str).eq("Gene Expression")
        if ge_mask.sum() > 0:
            adata = adata[:, ge_mask].copy()

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

    n_in_tissue = adata.n_obs

    adata = adata[
        (adata.obs["total_counts"] >= MIN_TOTAL_COUNTS)
        & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["pct_counts_mt"] < MAX_PCT_MT),
        :
    ].copy()

    # Spatial coordinates for Scanpy-compatible embedding.
    adata.obsm["X_spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)

    if scale_file.exists():
        with gzip.open(scale_file, "rt") as handle:
            adata.uns["spatial_scalefactors"] = json.load(handle)

    print("  Normalising/log-transforming and scoring signatures...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    available = set(adata.var_names)
    for sig_name, genes in signature_sets.items():
        present = [g for g in genes if g in available]
        adata.uns[f"{sig_name}_genes_present"] = present

        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=sig_name)
        else:
            adata.obs[sig_name] = 0.0

        sig_rows.append({
            "section_id": section_id,
            "signature": sig_name,
            "n_genes_requested": len(genes),
            "n_genes_present": len(present),
            "genes_present": ";".join(present),
            "mean_score": round(float(adata.obs[sig_name].mean()), 6) if adata.n_obs else np.nan,
            "max_score": round(float(adata.obs[sig_name].max()), 6) if adata.n_obs else np.nan,
        })

    out_file = OUT_DIR / f"{prefix}_spatial_scored.h5ad"
    adata.write_h5ad(out_file)

    summary_rows.append({
        "prefix": prefix,
        "section_id": section_id,
        "raw_spots": len(barcodes),
        "in_tissue_spots": n_in_tissue,
        "filtered_spots": adata.n_obs,
        "median_total_counts": round(float(adata.obs["total_counts"].median()), 2) if adata.n_obs else np.nan,
        "median_n_genes": round(float(adata.obs["n_genes_by_counts"].median()), 2) if adata.n_obs else np.nan,
        "median_pct_mt": round(float(adata.obs["pct_counts_mt"].median()), 3) if adata.n_obs else np.nan,
        "output_h5ad": str(out_file),
    })

    print(f"  Wrote {out_file} with {adata.n_obs} filtered tissue spots")

    # Basic first-pass spatial plots.
    for sig in [
        "sig_TEC",
        "sig_cTEC_like",
        "sig_mTEC_like",
        "sig_MHCII_AP",
        "sig_B_cell_APC",
        "sig_macrophage_APC",
        "sig_pDC_APC",
        "sig_fibroblast",
        "sig_endothelial",
    ]:
        plot_spatial_score(
            adata,
            section_id,
            sig,
            FIG_DIR / f"{prefix}_{sig}_spatial_score.png"
        )

    del adata, mat
    gc.collect()

summary = pd.DataFrame(summary_rows)
sig_summary = pd.DataFrame(sig_rows)

summary_out = META_DIR / "heimli_spatial_qc_summary.tsv"
sig_out = TAB_DIR / "heimli_spatial_signature_score_summary.tsv"

summary.to_csv(summary_out, sep="\t", index=False)
sig_summary.to_csv(sig_out, sep="\t", index=False)

print("Done.")
print(f"Wrote spatial QC summary: {summary_out}")
print(f"Wrote spatial signature summary: {sig_out}")
print()
print("Spatial QC summary:")
print(summary.to_string(index=False))
print()
print("Spatial signature score summary, first 40 rows:")
print(sig_summary.head(40).to_string(index=False))
