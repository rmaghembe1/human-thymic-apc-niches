#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import scanpy as sc

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")

IN_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_umap.h5ad"
ANNOT = PROJECT_DIR / "metadata/heimli_downsampled_cluster_annotations.tsv"
OUT_H5AD = PROJECT_DIR / "data/processed/heimli_gse207206_downsampled_annotated.h5ad"

FIG_DIR = PROJECT_DIR / "results/figures"
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

print("Reading Heimli downsampled UMAP object...")
adata = sc.read_h5ad(IN_H5AD)

print("Reading Heimli annotation table...")
annot = pd.read_csv(ANNOT, sep="\t")
annot["leiden_0_8"] = annot["leiden_0_8"].astype(str)

for col in ["sample_prefix", "donor", "fraction", "leiden_0_8"]:
    adata.obs[col] = adata.obs[col].astype(str)

mapping = annot.set_index("leiden_0_8")

adata.obs["working_label"] = adata.obs["leiden_0_8"].map(mapping["working_label"])
adata.obs["broad_class"] = adata.obs["leiden_0_8"].map(mapping["broad_class"])
adata.obs["annotation_confidence"] = adata.obs["leiden_0_8"].map(mapping["confidence"])

missing = adata.obs["working_label"].isna().sum()
if missing > 0:
    missing_clusters = sorted(adata.obs.loc[adata.obs["working_label"].isna(), "leiden_0_8"].unique())
    raise ValueError(f"Missing annotations for {missing} cells. Missing clusters: {missing_clusters}")

print("Writing annotated Heimli object...")
adata.write_h5ad(OUT_H5AD)

print("Building composition tables...")

overall_class = (
    adata.obs
    .groupby("broad_class", observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
overall_class["pct_total"] = (100 * overall_class["n_cells"] / adata.n_obs).round(2)

overall_cluster = (
    adata.obs
    .groupby(["leiden_0_8", "working_label", "broad_class", "annotation_confidence"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
overall_cluster["pct_total"] = (100 * overall_cluster["n_cells"] / adata.n_obs).round(2)

fraction_class = (
    adata.obs
    .groupby(["fraction", "broad_class"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
fraction_totals = (
    adata.obs
    .groupby("fraction", observed=True)
    .size()
    .rename("fraction_total_cells")
    .reset_index()
)
fraction_class = fraction_class.merge(fraction_totals, on="fraction", how="left")
fraction_class["pct_of_fraction"] = (
    100 * fraction_class["n_cells"] / fraction_class["fraction_total_cells"]
).round(2)

donor_class = (
    adata.obs
    .groupby(["donor", "fraction", "broad_class"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
donor_fraction_totals = (
    adata.obs
    .groupby(["donor", "fraction"], observed=True)
    .size()
    .rename("donor_fraction_total_cells")
    .reset_index()
)
donor_class = donor_class.merge(donor_fraction_totals, on=["donor", "fraction"], how="left")
donor_class["pct_of_donor_fraction"] = (
    100 * donor_class["n_cells"] / donor_class["donor_fraction_total_cells"]
).round(2)

overall_class_out = TAB_DIR / "heimli_downsampled_overall_broad_class_composition.tsv"
overall_cluster_out = TAB_DIR / "heimli_downsampled_overall_cluster_composition.tsv"
fraction_class_out = TAB_DIR / "heimli_downsampled_fraction_by_broad_class_composition.tsv"
donor_class_out = TAB_DIR / "heimli_downsampled_donor_fraction_by_broad_class_composition.tsv"

overall_class.to_csv(overall_class_out, sep="\t", index=False)
overall_cluster.to_csv(overall_cluster_out, sep="\t", index=False)
fraction_class.to_csv(fraction_class_out, sep="\t", index=False)
donor_class.to_csv(donor_class_out, sep="\t", index=False)

print("Plotting annotated UMAPs...")
sc.pl.umap(
    adata,
    color=["broad_class"],
    frameon=False,
    save="_heimli_annotated_broad_class.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["working_label"],
    frameon=False,
    save="_heimli_annotated_working_label.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["annotation_confidence"],
    frameon=False,
    save="_heimli_annotation_confidence.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["fraction", "donor"],
    frameon=False,
    save="_heimli_annotated_fraction_donor.png",
    show=False
)

print("Done.")
print(f"Wrote annotated object: {OUT_H5AD}")
print(f"Wrote: {overall_class_out}")
print(f"Wrote: {overall_cluster_out}")
print(f"Wrote: {fraction_class_out}")
print(f"Wrote: {donor_class_out}")
print()
print("Overall broad-class composition:")
print(overall_class.sort_values("n_cells", ascending=False).to_string(index=False))
print()
print("Fraction-by-broad-class composition:")
print(fraction_class.sort_values(["fraction", "broad_class"]).to_string(index=False))
