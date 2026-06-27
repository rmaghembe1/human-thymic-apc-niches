#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"
DOC_DIR.mkdir(parents=True, exist_ok=True)

in_md = DOC_DIR / "manuscript_draft_v1.md"
out_md = DOC_DIR / "manuscript_draft_v2_cited.md"
ref_tsv = DOC_DIR / "core_references_v1.tsv"

if not in_md.exists():
    raise FileNotFoundError(f"Missing manuscript draft: {in_md}")

refs = [
    {
        "key": "Bautista2021",
        "citation": "Bautista JL, Cramer NT, Miller CN, Chavez J, Berrios DI, Byrnes LE, et al. Single-cell transcriptional profiling of human thymic stroma uncovers novel cellular heterogeneity in the thymic medulla. Nature Communications. 2021;12:1096. doi:10.1038/s41467-021-21346-6.",
        "use_for": "Bautista stromal thymus dataset; human thymic stromal/TEC heterogeneity; fetal-postnatal-adult stromal scRNA-seq.",
        "status": "core_dataset_reference"
    },
    {
        "key": "GSE147520",
        "citation": "Gene Expression Omnibus accession GSE147520. Single-cell transcriptional profiling of human thymic stroma.",
        "use_for": "Raw Bautista GEO accession and data availability.",
        "status": "dataset_accession"
    },
    {
        "key": "Heimli2023",
        "citation": "Heimli M, et al. Multimodal human thymic profiling reveals trajectories and cellular milieu for T agonist selection. Frontiers in Immunology. 2023;13:1092028. doi:10.3389/fimmu.2022.1092028.",
        "use_for": "Heimli paediatric thymus CITE-seq and Visium spatial transcriptomics dataset.",
        "status": "core_dataset_reference"
    },
    {
        "key": "GSE207206",
        "citation": "Gene Expression Omnibus accession GSE207206. Multimodal human thymic profiling reveals trajectories and cellular milieu for T agonist selection.",
        "use_for": "Raw Heimli GEO accession; CITE-seq and spatial transcriptomics files.",
        "status": "dataset_accession"
    },
    {
        "key": "Park2020",
        "citation": "Park JE, Botting RA, Dominguez Conde C, Popescu DM, Lavaert M, Kunz DJ, et al. A cell atlas of human thymic development defines T cell repertoire formation. Science. 2020;367:eaay3224. doi:10.1126/science.aay3224.",
        "use_for": "Human thymic development single-cell context.",
        "status": "background_reference"
    },
    {
        "key": "Yayon2024",
        "citation": "Yayon N, et al. A spatial human thymus cell atlas mapped to a continuous tissue axis. Nature. 2024.",
        "use_for": "Recent human thymus spatial atlas context.",
        "status": "background_reference_to_verify_details"
    },
    {
        "key": "Klein2014",
        "citation": "Klein L, Kyewski B, Allen PM, Hogquist KA. Positive and negative selection of the T cell repertoire: what thymocytes see and do not see. Nature Reviews Immunology. 2014;14:377-391. doi:10.1038/nri3667.",
        "use_for": "Classical thymic positive and negative selection background.",
        "status": "background_reference"
    },
    {
        "key": "Giraud2024",
        "citation": "Giraud M, et al. Antigen presentation for central tolerance induction. Nature Reviews Immunology. 2024.",
        "use_for": "mTEC, thymic dendritic cell and thymic B-cell antigen presentation in central tolerance.",
        "status": "background_reference_to_verify_details"
    },
    {
        "key": "Yamano2015",
        "citation": "Yamano T, Steinert M, Klein L. Thymic B cells and central T cell tolerance. Frontiers in Immunology. 2015;6:376. doi:10.3389/fimmu.2015.00376.",
        "use_for": "Thymic B-cell biology and central tolerance background.",
        "status": "background_reference"
    },
    {
        "key": "Perera2015",
        "citation": "Perera J, et al. Thymic B cells are licensed to present self antigens for central T cell tolerance. Immunity. 2015.",
        "use_for": "Thymic B-cell antigen-presentation function.",
        "status": "background_reference_to_verify_details"
    },
    {
        "key": "Qi2022",
        "citation": "Qi Y, Zhang R, Lu Y, Zou X, Yang W. Aire and Fezf2, two regulators in medullary thymic epithelial cells, control autoimmune diseases by regulating tissue-specific antigens. Frontiers in Immunology. 2022;13:948259. doi:10.3389/fimmu.2022.948259.",
        "use_for": "AIRE/FEZF2 and mTEC tissue-restricted antigen expression background.",
        "status": "background_reference"
    },
    {
        "key": "TenXVisium2025",
        "citation": "10x Genomics. What is the spatial resolution and configuration of the capture area of the Visium v1 Gene Expression Slide? 10x Genomics Knowledge Base.",
        "use_for": "Visium spot size and capture-area configuration.",
        "status": "methods_background_reference"
    },
    {
        "key": "SpatialTechnologies2023",
        "citation": "Spatial transcriptomic technologies review. Use for general statement that Visium spot-level data are not single-cell resolution.",
        "use_for": "General spatial transcriptomics limitation and interpretation.",
        "status": "background_reference_to_verify_details"
    },
    {
        "key": "Scanpy2018",
        "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15. doi:10.1186/s13059-017-1382-0.",
        "use_for": "Scanpy software citation.",
        "status": "software_reference"
    },
    {
        "key": "UMAP2018",
        "citation": "McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426. 2018.",
        "use_for": "UMAP dimensionality reduction citation.",
        "status": "software_reference"
    },
    {
        "key": "Leiden2019",
        "citation": "Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports. 2019;9:5233. doi:10.1038/s41598-019-41695-z.",
        "use_for": "Leiden clustering citation.",
        "status": "software_reference"
    },
    {
        "key": "AnnData2024",
        "citation": "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. Journal of Open Source Software. 2024;9(101):4371. doi:10.21105/joss.04371.",
        "use_for": "AnnData software citation.",
        "status": "software_reference"
    },
]

pd.DataFrame(refs).to_csv(ref_tsv, sep="\t", index=False)

text = in_md.read_text()

replacements = {
    "The thymus is the central organ for T-cell development, selection and immune repertoire formation.":
        "The thymus is the central organ for T-cell development, selection and immune repertoire formation [@Klein2014].",

    "Classical models emphasise cortical thymic epithelial cells as key mediators of positive selection and medullary thymic epithelial cells as central organisers of negative selection and tissue-restricted antigen presentation.":
        "Classical models emphasise cortical thymic epithelial cells as key mediators of positive selection and medullary thymic epithelial cells as central organisers of negative selection and tissue-restricted antigen presentation [@Klein2014; @Qi2022].",

    "Dendritic cells, macrophages, B cells and other immune antigen-presenting populations also occupy the thymic environment and may contribute to the local organisation of tolerance and T-cell education.":
        "Dendritic cells, macrophages, B cells and other immune antigen-presenting populations also occupy the thymic environment and may contribute to the local organisation of tolerance and T-cell education [@Giraud2024; @Yamano2015; @Perera2015].",

    "Single-cell RNA sequencing has substantially advanced the resolution of human thymic cell-state maps.":
        "Single-cell RNA sequencing has substantially advanced the resolution of human thymic cell-state maps [@Park2020; @Bautista2021].",

    "Spatial transcriptomic methods such as Visium provide complementary information by projecting gene-expression programmes back into tissue sections, although spot-level resolution requires careful interpretation because each spot may contain multiple cells.":
        "Spatial transcriptomic methods such as Visium provide complementary information by projecting gene-expression programmes back into tissue sections, although spot-level resolution requires careful interpretation because each spot may contain multiple cells [@TenXVisium2025; @SpatialTechnologies2023].",

    "The Bautista human thymic stromal single-cell dataset was used as a reference for thymic epithelial and stromal architecture across fetal, postnatal and adult samples.":
        "The Bautista human thymic stromal single-cell dataset was used as a reference for thymic epithelial and stromal architecture across fetal, postnatal and adult samples [@Bautista2021; @GSE147520].",

    "The Heimli paediatric thymus CITE-seq and Visium dataset was used to resolve immune antigen-presenting, thymic B-cell, epithelial, stromal and spatially projected niche programmes.":
        "The Heimli paediatric thymus CITE-seq and Visium dataset was used to resolve immune antigen-presenting, thymic B-cell, epithelial, stromal and spatially projected niche programmes [@Heimli2023; @GSE207206].",

    "The Bautista dataset, GSE147520, contains human thymic stromal single-cell transcriptomes across fetal, postnatal and adult samples.":
        "The Bautista dataset, GSE147520, contains human thymic stromal single-cell transcriptomes across fetal, postnatal and adult samples [@Bautista2021; @GSE147520].",

    "The Heimli dataset, GSE207206, contains paediatric thymus CITE-seq and Visium spatial transcriptomic data.":
        "The Heimli dataset, GSE207206, contains paediatric thymus CITE-seq and Visium spatial transcriptomic data [@Heimli2023; @GSE207206].",

    "Analyses were performed in Python using Scanpy, AnnData, pandas, NumPy, SciPy and Matplotlib in a WSL/Linux environment.":
        "Analyses were performed in Python using Scanpy, AnnData, pandas, NumPy, SciPy and Matplotlib in a WSL/Linux environment [@Scanpy2018; @AnnData2024; @UMAP2018; @Leiden2019].",
}

for old, new in replacements.items():
    if old in text:
        text = text.replace(old, new)
    else:
        print(f"WARNING: phrase not found for replacement:\n{old}\n")

# Replace the placeholder References section with a core working list.
ref_lines = ["## References", ""]
for i, ref in enumerate(refs, start=1):
    ref_lines.append(f"{i}. [{ref['key']}] {ref['citation']}")
ref_block = "\n".join(ref_lines) + "\n"

if "## References\n\nReferences will be added and formatted according to the target journal style.\n" in text:
    text = text.replace(
        "## References\n\nReferences will be added and formatted according to the target journal style.\n",
        ref_block
    )
else:
    text = text.rstrip() + "\n\n" + ref_block

out_md.write_text(text)

print("Done.")
print(f"Wrote cited manuscript draft: {out_md}")
print(f"Wrote core reference table: {ref_tsv}")
print()
print("Next files:")
print(out_md)
print(ref_tsv)
