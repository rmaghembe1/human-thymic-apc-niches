#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import scanpy as sc

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")

IN_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_refined_scores.h5ad"
ANNOT = PROJECT_DIR / "metadata/bautista_downsampled_cluster_annotations.tsv"
OUT_H5AD = PROJECT_DIR / "data/processed/bautista_gse147520_downsampled_annotated.h5ad"

FIG_DIR = PROJECT_DIR / "results/figures"
TAB_DIR = PROJECT_DIR / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

print("Reading refined scored object...")
adata = sc.read_h5ad(IN_H5AD)

print("Reading cluster annotations...")
annot = pd.read_csv(ANNOT, sep="\t")
annot["leiden_0_6"] = annot["leiden_0_6"].astype(str)

for col in ["sample_prefix", "sample_short", "age_group", "age_value", "leiden_0_6"]:
    adata.obs[col] = adata.obs[col].astype(str)

mapping = annot.set_index("leiden_0_6")

adata.obs["working_label"] = adata.obs["leiden_0_6"].map(mapping["working_label"])
adata.obs["broad_class"] = adata.obs["leiden_0_6"].map(mapping["broad_class"])
adata.obs["annotation_confidence"] = adata.obs["leiden_0_6"].map(mapping["confidence"])

missing = adata.obs["working_label"].isna().sum()
if missing > 0:
    raise ValueError(f"{missing} cells have missing annotation labels. Check cluster annotation table.")

print("Writing annotated h5ad...")
adata.write_h5ad(OUT_H5AD)

print("Building composition tables...")

sample_class_counts = (
    adata.obs
    .groupby(["sample_short", "age_group", "age_value", "broad_class"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)

sample_totals = (
    adata.obs
    .groupby(["sample_short"], observed=True)
    .size()
    .rename("sample_total_cells")
    .reset_index()
)

sample_class_counts = sample_class_counts.merge(sample_totals, on="sample_short", how="left")
sample_class_counts["pct_of_sample"] = (
    100 * sample_class_counts["n_cells"] / sample_class_counts["sample_total_cells"]
).round(2)

sample_cluster_counts = (
    adata.obs
    .groupby(["sample_short", "age_group", "age_value", "leiden_0_6", "working_label", "broad_class"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)

sample_cluster_counts = sample_cluster_counts.merge(sample_totals, on="sample_short", how="left")
sample_cluster_counts["pct_of_sample"] = (
    100 * sample_cluster_counts["n_cells"] / sample_cluster_counts["sample_total_cells"]
).round(2)

overall_class_counts = (
    adata.obs
    .groupby(["broad_class"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
overall_class_counts["pct_total"] = (
    100 * overall_class_counts["n_cells"] / adata.n_obs
).round(2)

overall_cluster_counts = (
    adata.obs
    .groupby(["leiden_0_6", "working_label", "broad_class", "annotation_confidence"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
overall_cluster_counts["pct_total"] = (
    100 * overall_cluster_counts["n_cells"] / adata.n_obs
).round(2)

sample_class_out = TAB_DIR / "bautista_downsampled_sample_by_broad_class_composition.tsv"
sample_cluster_out = TAB_DIR / "bautista_downsampled_sample_by_cluster_composition.tsv"
overall_class_out = TAB_DIR / "bautista_downsampled_overall_broad_class_composition.tsv"
overall_cluster_out = TAB_DIR / "bautista_downsampled_overall_cluster_composition.tsv"

sample_class_counts.to_csv(sample_class_out, sep="\t", index=False)
sample_cluster_counts.to_csv(sample_cluster_out, sep="\t", index=False)
overall_class_counts.to_csv(overall_class_out, sep="\t", index=False)
overall_cluster_counts.to_csv(overall_cluster_out, sep="\t", index=False)

print("Plotting annotated UMAPs...")
sc.pl.umap(
    adata,
    color=["broad_class"],
    frameon=False,
    save="_bautista_annotated_broad_class.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["working_label"],
    frameon=False,
    save="_bautista_annotated_working_label.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["annotation_confidence"],
    frameon=False,
    save="_bautista_annotation_confidence.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["sample_short", "age_group"],
    frameon=False,
    save="_bautista_annotated_sample_age.png",
    show=False
)

print("Done.")
print(f"Wrote annotated object: {OUT_H5AD}")
print(f"Wrote: {sample_class_out}")
print(f"Wrote: {sample_cluster_out}")
print(f"Wrote: {overall_class_out}")
print(f"Wrote: {overall_cluster_out}")
print()
print("Overall broad-class composition:")
print(overall_class_counts.sort_values("n_cells", ascending=False).to_string(index=False))
print()
print("Sample-by-broad-class composition:")
print(sample_class_counts.sort_values(["sample_short", "broad_class"]).to_string(index=False))
