#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"
FIG_DIR = PROJECT_DIR / "results/figures/manuscript_final_panels"
TAB_DIR = PROJECT_DIR / "results/tables"
META_DIR = PROJECT_DIR / "metadata"
SCRIPT_DIR = PROJECT_DIR / "scripts"

DOC_DIR.mkdir(parents=True, exist_ok=True)

items = []

def add_item(category, label, path, status_note):
    path = Path(path)
    items.append({
        "category": category,
        "label": label,
        "path": str(path),
        "exists": path.exists(),
        "size_kb": round(path.stat().st_size / 1024, 1) if path.exists() else None,
        "status_note": status_note,
    })

# Main manuscript files
add_item(
    "manuscript",
    "Reference-checked manuscript draft",
    DOC_DIR / "manuscript_draft_v3_reference_checked.md",
    "Current working manuscript draft"
)

add_item(
    "references",
    "Checked core reference table",
    DOC_DIR / "core_references_v2_checked.tsv",
    "Reference table with verification status"
)

add_item(
    "references",
    "Reference verification notes",
    DOC_DIR / "reference_verification_notes_v2.md",
    "Items still requiring final verification before submission"
)

# Main figures
for fig in ["figure1_final_dataset_architecture", "figure2_final_marker_programmes", "figure3_final_spatial_copatterns"]:
    for ext in ["png", "svg", "pdf"]:
        add_item(
            "main_figure",
            f"{fig}.{ext}",
            FIG_DIR / f"{fig}.{ext}",
            "Frozen or near-final main figure asset"
        )

# Core supplementary / results tables
core_tables = [
    ("Bautista QC summary", META_DIR / "bautista_gse147520_qc_summary_clean.tsv"),
    ("Bautista cluster annotations", META_DIR / "bautista_downsampled_cluster_annotations.tsv"),
    ("Bautista broad-class composition", TAB_DIR / "bautista_downsampled_overall_broad_class_composition.tsv"),
    ("Bautista cluster composition", TAB_DIR / "bautista_downsampled_overall_cluster_composition.tsv"),
    ("Bautista marker score summary", TAB_DIR / "bautista_downsampled_refined_marker_score_summary.tsv"),
    ("Bautista RNA markers", TAB_DIR / "bautista_downsampled_rank_genes_by_leiden.tsv"),
    ("Heimli per-sample QC", META_DIR / "heimli_gse207206_filtered_per_sample_qc_summary.tsv"),
    ("Heimli cluster annotations", META_DIR / "heimli_downsampled_cluster_annotations.tsv"),
    ("Heimli broad-class composition", TAB_DIR / "heimli_downsampled_overall_broad_class_composition.tsv"),
    ("Heimli fraction composition", TAB_DIR / "heimli_downsampled_fraction_by_broad_class_composition.tsv"),
    ("Heimli RNA-only markers", TAB_DIR / "heimli_downsampled_rank_genes_by_leiden_RNA_only.tsv"),
    ("Heimli ADT mean by cluster", TAB_DIR / "heimli_downsampled_ADT_mean_by_leiden.tsv"),
    ("Integrated evidence summary", TAB_DIR / "integrated_bautista_heimli_evidence_summary.tsv"),
    ("Heimli spatial QC", META_DIR / "heimli_spatial_qc_summary.tsv"),
    ("Heimli spatial signature summary", TAB_DIR / "heimli_spatial_signature_score_summary.tsv"),
    ("Heimli spatial ranked signatures", TAB_DIR / "heimli_spatial_signature_ranked_by_section.tsv"),
    ("Heimli spatial co-high summary", TAB_DIR / "heimli_spatial_cohigh_summary_across_sections.tsv"),
    ("Heimli spatial focus correlations", TAB_DIR / "heimli_spatial_focus_correlation_summary_across_sections.tsv"),
]

for label, path in core_tables:
    add_item(
        "supplementary_table_candidate",
        label,
        path,
        "Candidate supplementary table or source table"
    )

# Core scripts
core_scripts = [
    "build_bautista_filtered_h5ad.py",
    "process_bautista_downsampled_umap.py",
    "refine_bautista_downsampled_marker_scores.py",
    "bautista_downsampled_rank_genes.py",
    "apply_bautista_annotations_and_composition.py",
    "build_heimli_filtered_per_sample_h5ad.py",
    "process_heimli_downsampled_umap.py",
    "heimli_downsampled_rank_RNA_and_ADT_separate.py",
    "apply_heimli_annotations_and_composition.py",
    "process_heimli_spatial_signature_scores.py",
    "summarise_heimli_spatial_signatures.py",
    "analyse_heimli_spatial_copatterns.py",
    "write_heimli_spatial_copattern_summary.py",
    "make_figure1_final_layout.py",
    "make_figure2_final_layout.py",
    "make_figure3_final_layout.py",
    "write_manuscript_draft_v1.py",
    "add_core_citations_to_manuscript.py",
    "check_and_clean_references_v3.py",
]

for script in core_scripts:
    add_item(
        "analysis_script",
        script,
        SCRIPT_DIR / script,
        "Core reproducibility script"
    )

manifest = pd.DataFrame(items)
manifest_out = DOC_DIR / "submission_readiness_manifest_v1.tsv"
manifest.to_csv(manifest_out, sep="\t", index=False)

missing = manifest[manifest["exists"] == False].copy()
present = manifest[manifest["exists"] == True].copy()

# Supplementary material plan
supp_lines = []
supp_lines.append("# Supplementary material plan v1\n")
supp_lines.append("## Supplementary tables\n")
supp_lines.append("| Proposed item | Source file | Purpose |")
supp_lines.append("|---|---|---|")

supp_table_map = [
    ("Table S1", "bautista_gse147520_qc_summary_clean.tsv", "Bautista sample-level QC summary."),
    ("Table S2", "bautista_downsampled_cluster_annotations.tsv", "Bautista exploratory cluster annotation table."),
    ("Table S3", "bautista_downsampled_refined_marker_score_summary.tsv", "Bautista refined marker-score summary by cluster."),
    ("Table S4", "bautista_downsampled_rank_genes_by_leiden.tsv", "Bautista ranked marker genes by cluster."),
    ("Table S5", "heimli_gse207206_filtered_per_sample_qc_summary.tsv", "Heimli per-library QC summary."),
    ("Table S6", "heimli_downsampled_cluster_annotations.tsv", "Heimli exploratory cluster annotation table."),
    ("Table S7", "heimli_downsampled_rank_genes_by_leiden_RNA_only.tsv", "Heimli RNA-only ranked marker genes by cluster."),
    ("Table S8", "heimli_downsampled_ADT_mean_by_leiden.tsv", "Heimli ADT/protein mean values by cluster."),
    ("Table S9", "integrated_bautista_heimli_evidence_summary.tsv", "Integrated evidence-axis summary and claim boundaries."),
    ("Table S10", "heimli_spatial_qc_summary.tsv", "Heimli Visium section-level QC summary."),
    ("Table S11", "heimli_spatial_signature_ranked_by_section.tsv", "Spatial signature ranking by section."),
    ("Table S12", "heimli_spatial_cohigh_summary_across_sections.tsv", "Spatial co-high summary across sections."),
    ("Table S13", "heimli_spatial_focus_correlation_summary_across_sections.tsv", "Spatial focus-correlation summary across sections."),
]

for item, source, purpose in supp_table_map:
    supp_lines.append(f"| {item} | `{source}` | {purpose} |")

supp_lines.append("\n## Supplementary figures\n")
supp_lines.append("| Proposed item | Purpose |")
supp_lines.append("|---|---|")
supp_lines.append("| Figure S1 | Bautista detailed working-label UMAP and full refined marker-score grid. |")
supp_lines.append("| Figure S2 | Heimli full detailed broad-class UMAP and complete marker-score grid. |")
supp_lines.append("| Figure S3 | Cross-section spatial maps for MHC-II/APC, B-cell APC, mTEC-like, cTEC-like and thymocyte signatures. |")
supp_lines.append("| Figure S4 | Spatial dominant-axis maps and full co-pattern support panels. |")

supp_out = DOC_DIR / "supplementary_material_plan_v1.md"
supp_out.write_text("\n".join(supp_lines) + "\n")

# Readiness report
report_lines = []
report_lines.append("# Submission readiness report v1\n")
report_lines.append("## Current status\n")
report_lines.append("- Main manuscript draft: `docs/manuscript_draft_v3_reference_checked.md`.")
report_lines.append("- Main figures 1-3: frozen / near-final in PNG, SVG and PDF formats.")
report_lines.append("- Core reference table: `docs/core_references_v2_checked.tsv`.")
report_lines.append("- Supplementary material plan: `docs/supplementary_material_plan_v1.md`.")
report_lines.append("\n## File completeness\n")
report_lines.append(f"- Present files in manifest: {len(present)}")
report_lines.append(f"- Missing files in manifest: {len(missing)}")

if len(missing) > 0:
    report_lines.append("\n## Missing or not yet generated files\n")
    report_lines.append(missing[["category", "label", "path"]].to_markdown(index=False))
else:
    report_lines.append("\nNo missing files were detected among tracked core items.")

report_lines.append("\n## Recommended next actions\n")
report_lines.append("1. Verify remaining reference items marked as accession or keep-verify-before-submission.")
report_lines.append("2. Generate Supplementary Figures S1-S4 from existing outputs.")
report_lines.append("3. Convert the manuscript from Markdown to DOCX/PDF only after reference and supplementary material checks.")
report_lines.append("4. Decide target journal and adjust format, word count, section order and citation style.")
report_lines.append("5. Prepare data/code availability with a repository structure and README.")

report_out = DOC_DIR / "submission_readiness_report_v1.md"
report_out.write_text("\n".join(report_lines) + "\n")

print("Done.")
print(f"Wrote manifest: {manifest_out}")
print(f"Wrote supplementary material plan: {supp_out}")
print(f"Wrote readiness report: {report_out}")
print()
print("Present files:", len(present))
print("Missing files:", len(missing))
if len(missing) > 0:
    print()
    print("Missing tracked files:")
    print(missing[["category", "label", "path"]].to_string(index=False))
