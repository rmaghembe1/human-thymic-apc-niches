#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
FIG_DIR = PROJECT_DIR / "results/figures"
SPATIAL_COPATTERN_DIR = FIG_DIR / "heimli_spatial/copatterns"
OUT_DIR = FIG_DIR / "supplementary_figures"
DOC_DIR = PROJECT_DIR / "docs"

OUT_DIR.mkdir(parents=True, exist_ok=True)
DOC_DIR.mkdir(parents=True, exist_ok=True)

def read_img(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing image: {path}")
    return mpimg.imread(path)

def save_fig_all(fig, basename):
    png = OUT_DIR / f"{basename}.png"
    svg = OUT_DIR / f"{basename}.svg"
    pdf = OUT_DIR / f"{basename}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png}")
    print(f"Wrote {svg}")
    print(f"Wrote {pdf}")

def make_contact_sheet(panels, title, basename, ncols=2, figsize_per_panel=(6.0, 4.8)):
    existing = []
    missing = []
    for label, path in panels:
        path = Path(path)
        if path.exists():
            existing.append((label, path))
        else:
            missing.append((label, path))

    if not existing:
        raise FileNotFoundError(f"No panels found for {basename}")

    nrows = (len(existing) + ncols - 1) // ncols
    fig = plt.figure(figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))

    for i, (label, path) in enumerate(existing, start=1):
        ax = fig.add_subplot(nrows, ncols, i)
        ax.imshow(read_img(path))
        ax.axis("off")
        ax.set_title(label, fontsize=11, loc="left")

    fig.suptitle(title, fontsize=15, y=0.995)
    plt.tight_layout()
    save_fig_all(fig, basename)

    return missing

all_missing = []

# ----------------------------
# Supplementary Figure S1
# Bautista detailed stromal/TEC annotation and marker evidence
# ----------------------------
s1_panels = [
    ("A. Bautista broad-class annotation", FIG_DIR / "umap_bautista_annotated_broad_class.png"),
    ("B. Bautista working-label annotation", FIG_DIR / "umap_bautista_annotated_working_label.png"),
    ("C. Bautista annotation confidence", FIG_DIR / "umap_bautista_annotation_confidence.png"),
    ("D. Bautista refined marker-score projections", FIG_DIR / "umap_bautista_downsampled_refined_scores.png"),
]

missing = make_contact_sheet(
    s1_panels,
    "Supplementary Figure S1. Bautista stromal/TEC annotation and marker evidence",
    "figureS1_bautista_stromal_TEC_annotation_marker_evidence",
    ncols=2,
    figsize_per_panel=(6.2, 5.0)
)
all_missing.extend([("Figure S1", x[0], str(x[1])) for x in missing])

# ----------------------------
# Supplementary Figure S2
# Heimli detailed annotation and marker evidence
# ----------------------------
s2_panels = [
    ("A. Heimli broad-class annotation", FIG_DIR / "umap_heimli_annotated_broad_class.png"),
    ("B. Heimli working-label annotation", FIG_DIR / "umap_heimli_annotated_working_label.png"),
    ("C. Heimli annotation confidence", FIG_DIR / "umap_heimli_annotation_confidence.png"),
    ("D. Heimli full marker-score projections", FIG_DIR / "umap_heimli_downsampled_marker_scores.png"),
    ("E. Heimli fraction/donor projection", FIG_DIR / "umap_heimli_annotated_fraction_donor.png"),
]

missing = make_contact_sheet(
    s2_panels,
    "Supplementary Figure S2. Heimli detailed annotation and marker evidence",
    "figureS2_heimli_detailed_annotation_marker_evidence",
    ncols=2,
    figsize_per_panel=(6.2, 5.0)
)
all_missing.extend([("Figure S2", x[0], str(x[1])) for x in missing])

# ----------------------------
# Supplementary Figure S3
# Cross-section spatial maps for key signatures
# ----------------------------
sections = ["S1_A1", "S1_B1", "S1_C1", "S1_D1", "S2_A1", "S2_B1", "S2_C1", "S2_D1"]

sig_defs = [
    ("MHC-II/APC", "sig_MHCII_AP_z"),
    ("B-cell APC", "sig_B_cell_APC_z"),
    ("mTEC-like", "sig_mTEC_like_z"),
    ("cTEC-like", "sig_cTEC_like_z"),
    ("T-cell/thymocyte", "axis_T_cell_thymocyte_z"),
]

for sig_label, sig_key in sig_defs:
    panels = []
    for section in sections:
        panels.append((section, SPATIAL_COPATTERN_DIR / f"{section}_{sig_key}.png"))

    safe_label = sig_label.replace("/", "_").replace(" ", "_").replace("-", "_").lower()
    missing = make_contact_sheet(
        panels,
        f"Supplementary Figure S3. Cross-section spatial maps: {sig_label}",
        f"figureS3_spatial_cross_section_{safe_label}",
        ncols=4,
        figsize_per_panel=(3.6, 3.6)
    )
    all_missing.extend([(f"Figure S3 {sig_label}", x[0], str(x[1])) for x in missing])

# ----------------------------
# Supplementary Figure S4
# Dominant-axis and dominant-signature spatial support
# ----------------------------
s4_panels = []
for section in sections:
    s4_panels.append((f"{section} dominant axis", SPATIAL_COPATTERN_DIR / f"{section}_dominant_axis.png"))

missing = make_contact_sheet(
    s4_panels,
    "Supplementary Figure S4A. Spatial dominant-axis maps across Heimli sections",
    "figureS4A_spatial_dominant_axis_maps",
    ncols=4,
    figsize_per_panel=(3.6, 3.6)
)
all_missing.extend([("Figure S4A", x[0], str(x[1])) for x in missing])

s4b_panels = []
for section in sections:
    s4b_panels.append((f"{section} dominant signature", SPATIAL_COPATTERN_DIR / f"{section}_dominant_signature.png"))

missing = make_contact_sheet(
    s4b_panels,
    "Supplementary Figure S4B. Spatial dominant-signature maps across Heimli sections",
    "figureS4B_spatial_dominant_signature_maps",
    ncols=4,
    figsize_per_panel=(3.6, 3.6)
)
all_missing.extend([("Figure S4B", x[0], str(x[1])) for x in missing])

# ----------------------------
# Write supplementary legend draft
# ----------------------------
legend_text = """# Supplementary figure legends v1

**Supplementary Figure S1. Bautista stromal/TEC annotation and marker evidence.**  
Supplementary UMAP projections showing Bautista broad-class annotation, working-label annotation, annotation confidence and refined marker-score projections. These panels support use of the Bautista dataset as a stromal/TEC architecture reference and document the exploratory annotation structure used for downstream interpretation.

**Supplementary Figure S2. Heimli detailed annotation and marker evidence.**  
Supplementary UMAP projections showing Heimli broad-class annotation, working-label annotation, annotation confidence, full marker-score projections and fraction/donor structure. These panels support the detailed annotation of immune APC, TEC-like, stromal, endothelial, thymocyte/T-cell and mixed low-confidence populations.

**Supplementary Figure S3. Cross-section spatial signature maps.**  
Cross-section spatial maps for key Heimli Visium signatures, including MHC-II/APC, B-cell APC, mTEC-like, cTEC-like and T-cell/thymocyte programmes. These maps support the reproducibility of the spatial patterns summarised in Figure 3.

**Supplementary Figure S4. Dominant spatial-axis and dominant-signature maps.**  
Spatial maps showing dominant niche-axis and dominant-signature assignments across Heimli Visium sections. These maps provide additional support for spatial partitioning of immune APC, epithelial, stromal and thymocyte-associated programmes.
"""

legend_out = DOC_DIR / "supplementary_figure_legends_v1.md"
legend_out.write_text(legend_text)

# ----------------------------
# Write missing-panel report
# ----------------------------
missing_out = DOC_DIR / "supplementary_figure_missing_panels_v1.tsv"

if all_missing:
    import pandas as pd
    pd.DataFrame(all_missing, columns=["supplementary_figure", "panel", "missing_path"]).to_csv(
        missing_out, sep="\t", index=False
    )
else:
    missing_out.write_text("No missing supplementary figure panels detected.\n")

print("Done.")
print(f"Wrote supplementary figures to: {OUT_DIR}")
print(f"Wrote supplementary legends: {legend_out}")
print(f"Wrote missing-panel report: {missing_out}")
print(f"Missing panels: {len(all_missing)}")
