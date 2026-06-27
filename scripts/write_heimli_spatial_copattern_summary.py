#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
TAB_DIR = PROJECT_DIR / "results/tables"
DOC_DIR = PROJECT_DIR / "docs"
DOC_DIR.mkdir(parents=True, exist_ok=True)

axis = pd.read_csv(TAB_DIR / "heimli_spatial_dominant_axis_by_section.tsv", sep="\t")
cohigh = pd.read_csv(TAB_DIR / "heimli_spatial_cohigh_rules_by_section.tsv", sep="\t")
corr = pd.read_csv(TAB_DIR / "heimli_spatial_focus_signature_correlations.tsv", sep="\t")

# Summaries across sections
axis_summary = (
    axis
    .groupby("dominant_axis", observed=True)
    .agg(
        mean_pct_spots=("pct_spots", "mean"),
        min_pct_spots=("pct_spots", "min"),
        max_pct_spots=("pct_spots", "max"),
        total_spots=("n_spots", "sum"),
    )
    .reset_index()
    .sort_values("mean_pct_spots", ascending=False)
)

cohigh_summary = (
    cohigh
    .groupby("cohigh_rule", observed=True)
    .agg(
        mean_pct_spots=("pct_spots", "mean"),
        min_pct_spots=("pct_spots", "min"),
        max_pct_spots=("pct_spots", "max"),
        total_cohigh_spots=("n_spots", "sum"),
    )
    .reset_index()
    .sort_values("mean_pct_spots", ascending=False)
)

corr_summary = (
    corr
    .assign(pair=lambda d: d["signature_1"] + " :: " + d["signature_2"])
    .groupby("pair", observed=True)
    .agg(
        mean_spearman_rho=("spearman_rho", "mean"),
        min_spearman_rho=("spearman_rho", "min"),
        max_spearman_rho=("spearman_rho", "max"),
        n_sections=("section_id", "nunique"),
    )
    .reset_index()
    .sort_values("mean_spearman_rho", ascending=False)
)

axis_out = TAB_DIR / "heimli_spatial_dominant_axis_summary_across_sections.tsv"
cohigh_out = TAB_DIR / "heimli_spatial_cohigh_summary_across_sections.tsv"
corr_out = TAB_DIR / "heimli_spatial_focus_correlation_summary_across_sections.tsv"

axis_summary.to_csv(axis_out, sep="\t", index=False)
cohigh_summary.to_csv(cohigh_out, sep="\t", index=False)
corr_summary.to_csv(corr_out, sep="\t", index=False)

summary_md = DOC_DIR / "heimli_spatial_copattern_interpretation.md"

lines = []
lines.append("# Heimli spatial co-pattern interpretation\n")
lines.append("## Core result\n")
lines.append("First-pass Visium signature co-pattern analysis identifies recurrent MHC-II/APC-rich spatial regions that co-enrich with B-cell/APC, mTEC-like, macrophage/APC and pDC-like signatures, while showing relative spatial separation from thymocyte and cTEC-like programmes.\n")

lines.append("## Dominant spatial axes across sections\n")
lines.append(axis_summary.to_markdown(index=False))
lines.append("\n")

lines.append("## Co-high spot patterns across sections\n")
lines.append(cohigh_summary.to_markdown(index=False))
lines.append("\n")

lines.append("## Focus spatial correlations across sections\n")
lines.append(corr_summary.to_markdown(index=False))
lines.append("\n")

lines.append("## Manuscript-ready interpretation\n")
lines.append("- B-cell/APC and MHC-II signatures show the most consistent positive spatial association across all sections.")
lines.append("- mTEC-like and MHC-II signatures also show recurrent positive association, supporting an epithelial antigen-presenting spatial component.")
lines.append("- Macrophage/APC and pDC-like signatures contribute additional immune APC-associated spatial programmes.")
lines.append("- Thymocyte and cTEC-like signatures are negatively correlated with the MHC-II/APC axis, supporting spatial partitioning of thymocyte/cortical-like and APC-rich regions.")
lines.append("- These results should be described as signature co-patterns, not direct cell-cell contact or formal deconvolution.")

lines.append("\n## Suggested Figure 3 structure\n")
lines.append("- Panel A: Heimli spatial section QC and workflow.")
lines.append("- Panel B: representative spatial maps for MHC-II/APC, B-cell/APC, mTEC-like, cTEC-like and thymocyte signatures.")
lines.append("- Panel C: dominant spatial axis maps.")
lines.append("- Panel D: co-high spot quantification across sections.")
lines.append("- Panel E: focus correlation heatmap showing MHC-II/APC relationships.")
lines.append("- Panel F: schematic model of spatially partitioned thymic antigen-presenting niches.")

summary_md.write_text("\n".join(lines))

print("Done.")
print(f"Wrote axis summary: {axis_out}")
print(f"Wrote co-high summary: {cohigh_out}")
print(f"Wrote correlation summary: {corr_out}")
print(f"Wrote interpretation document: {summary_md}")
print()
print("Co-high summary across sections:")
print(cohigh_summary.to_string(index=False))
print()
print("Focus correlation summary across sections:")
print(corr_summary.to_string(index=False))
