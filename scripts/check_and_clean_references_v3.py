#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"

in_md = DOC_DIR / "manuscript_draft_v2_cited.md"
out_md = DOC_DIR / "manuscript_draft_v3_reference_checked.md"
out_refs = DOC_DIR / "core_references_v2_checked.tsv"
notes_out = DOC_DIR / "reference_verification_notes_v2.md"

if not in_md.exists():
    raise FileNotFoundError(f"Missing manuscript: {in_md}")

refs = [
    {
        "key": "Bautista2021",
        "citation": "Bautista JL, Cramer NT, Miller CN, Chavez J, Berrios DI, Byrnes LE, Germino J, Ntranos V, Sneddon JB, Burt TD, Gardner JM, Ye CJ, Anderson MS, Parent AV. Single-cell transcriptional profiling of human thymic stroma uncovers novel cellular heterogeneity in the thymic medulla. Nature Communications. 2021;12:1096. doi:10.1038/s41467-021-21346-6.",
        "status": "verified",
        "note": "Primary Bautista human thymic stromal scRNA-seq dataset reference."
    },
    {
        "key": "GSE147520",
        "citation": "Gene Expression Omnibus accession GSE147520. Single-cell transcriptional profiling of human thymic stroma.",
        "status": "accession_to_confirm_before_submission",
        "note": "Raw Bautista GEO accession used for data retrieval."
    },
    {
        "key": "Heimli2023",
        "citation": "Heimli M, Flam ST, Hjorthaug HS, Trinh D, Frisk M, Dumont KA, Ribarska T, Tekpli X, Saare M, Lie BA. Multimodal human thymic profiling reveals trajectories and cellular milieu for T agonist selection. Frontiers in Immunology. 2023;13:1092028. doi:10.3389/fimmu.2022.1092028.",
        "status": "verified",
        "note": "Primary Heimli paediatric thymus CITE-seq and spatial transcriptomics dataset reference."
    },
    {
        "key": "GSE207206",
        "citation": "Gene Expression Omnibus accession GSE207206. Multimodal human thymic profiling reveals trajectories and cellular milieu for T agonist selection.",
        "status": "accession_to_confirm_before_submission",
        "note": "Raw Heimli GEO accession used for data retrieval."
    },
    {
        "key": "Park2020",
        "citation": "Park JE, Botting RA, Dominguez Conde C, Popescu DM, Lavaert M, Kunz DJ, et al. A cell atlas of human thymic development defines T cell repertoire formation. Science. 2020;367:eaay3224. doi:10.1126/science.aay3224.",
        "status": "verified",
        "note": "Human thymic development single-cell background."
    },
    {
        "key": "Yayon2024",
        "citation": "Yayon N, Kedlian VR, Boehme L, et al. A spatial human thymus cell atlas mapped to a continuous tissue axis. Nature. 2024;635:708-718. doi:10.1038/s41586-024-07944-6.",
        "status": "verified",
        "note": "Human thymus spatial atlas background."
    },
    {
        "key": "Klein2014",
        "citation": "Klein L, Kyewski B, Allen PM, Hogquist KA. Positive and negative selection of the T cell repertoire: what thymocytes see and do not see. Nature Reviews Immunology. 2014;14:377-391. doi:10.1038/nri3667.",
        "status": "verified",
        "note": "Classical thymic positive and negative selection background."
    },
    {
        "key": "KleinPetrozziello2025",
        "citation": "Klein L, Petrozziello E. Antigen presentation for central tolerance induction. Nature Reviews Immunology. 2025;25:57-72. doi:10.1038/s41577-024-01076-8.",
        "status": "verified",
        "note": "Current central tolerance and thymic APC review; replaces earlier provisional Giraud2024 key."
    },
    {
        "key": "Yamano2015Review",
        "citation": "Yamano T, Steinert M, Klein L. Thymic B cells and central T cell tolerance. Frontiers in Immunology. 2015;6:376. doi:10.3389/fimmu.2015.00376.",
        "status": "verified",
        "note": "Thymic B-cell central tolerance review."
    },
    {
        "key": "Yamano2015Immunity",
        "citation": "Yamano T, Nedjic J, Hinterberger M, Steinert M, Koser S, Pinto S, Gerdes N, Lutgens E, Ishimaru N, Busslinger M, Brors B, Kyewski B, Klein L. Thymic B cells are licensed to present self antigens for central T cell tolerance induction. Immunity. 2015;42:1048-1061. doi:10.1016/j.immuni.2015.05.013.",
        "status": "verified_from_secondary_and_source_lookup",
        "note": "Functional thymic B-cell antigen-presentation reference."
    },
    {
        "key": "Perera2013",
        "citation": "Perera J, Meng L, Meng F, Huang H. Autoreactive thymic B cells are efficient antigen-presenting cells of cognate self-antigens for T cell negative selection. Proceedings of the National Academy of Sciences of the United States of America. 2013;110:17011-17016.",
        "status": "verified_from_review_reference",
        "note": "Thymic B-cell negative selection reference."
    },
    {
        "key": "Perera2016",
        "citation": "Perera J, et al. Self-antigen-driven thymic B cell class switching promotes T cell central tolerance. Cell Reports. 2016;17:387-398.",
        "status": "verified_from_review_reference",
        "note": "Thymic B-cell class switching and central tolerance reference."
    },
    {
        "key": "Qi2022",
        "citation": "Qi Y, Zhang R, Lu Y, Zou X, Yang W. Aire and Fezf2, two regulators in medullary thymic epithelial cells, control autoimmune diseases by regulating tissue-specific antigens. Frontiers in Immunology. 2022;13:948259. doi:10.3389/fimmu.2022.948259.",
        "status": "keep_verify_before_submission",
        "note": "AIRE/FEZF2 and mTEC tissue-restricted antigen expression background."
    },
    {
        "key": "TenXVisium2025",
        "citation": "10x Genomics. What is the spatial resolution and configuration of the capture area of the Visium v1 Gene Expression Slide? 10x Genomics Knowledge Base.",
        "status": "verified",
        "note": "Visium v1 capture area, spot count, 55 micrometer spot diameter and 100 micrometer centre-to-centre spacing."
    },
    {
        "key": "Moffitt2022",
        "citation": "Moffitt JR, Lundberg E, Heyn H. The emerging landscape of spatial profiling technologies. Nature Reviews Genetics. 2022;23:741-759. doi:10.1038/s41576-022-00515-3.",
        "status": "verified_from_nature_reference_list",
        "note": "General spatial profiling technologies background; replaces vague SpatialTechnologies2023 placeholder."
    },
    {
        "key": "Scanpy2018",
        "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15. doi:10.1186/s13059-017-1382-0.",
        "status": "verified",
        "note": "Scanpy software citation."
    },
    {
        "key": "UMAP2018",
        "citation": "McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426. 2018.",
        "status": "standard_software_reference",
        "note": "UMAP dimensionality reduction citation."
    },
    {
        "key": "Leiden2019",
        "citation": "Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports. 2019;9:5233. doi:10.1038/s41598-019-41695-z.",
        "status": "verified",
        "note": "Leiden clustering citation."
    },
    {
        "key": "AnnData2024",
        "citation": "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. Journal of Open Source Software. 2024;9(101):4371. doi:10.21105/joss.04371.",
        "status": "keep_verify_before_submission",
        "note": "AnnData software citation."
    },
]

df = pd.DataFrame(refs)
df.to_csv(out_refs, sep="\t", index=False)

text = in_md.read_text()

# Replace obsolete provisional citation keys.
text = text.replace("@Giraud2024", "@KleinPetrozziello2025")
text = text.replace("@Perera2015", "@Yamano2015Immunity")
text = text.replace("@SpatialTechnologies2023", "@Moffitt2022")

# Replace old reference section with checked references.
ref_lines = ["## References", ""]
for i, row in enumerate(refs, start=1):
    ref_lines.append(f"{i}. [{row['key']}] {row['citation']}")
ref_block = "\n".join(ref_lines) + "\n"

if "## References\n" in text:
    text = text.split("## References\n")[0].rstrip() + "\n\n" + ref_block
else:
    text = text.rstrip() + "\n\n" + ref_block

out_md.write_text(text)

notes = """# Reference verification notes v2

## Checked and corrected

- Replaced provisional `Giraud2024` citation key with `KleinPetrozziello2025`.
- Replaced vague `SpatialTechnologies2023` placeholder with `Moffitt2022`.
- Replaced provisional `Perera2015` key with `Yamano2015Immunity`; retained Perera 2013 and Perera 2016 as supporting thymic B-cell tolerance references.
- Added full Yayon et al. Nature citation: Nature 635, 708-718 (2024), doi:10.1038/s41586-024-07944-6.
- Kept GEO accession entries as dataset-accession placeholders to confirm against GEO pages before submission.

## Still to verify before final submission

- GEO accession metadata formatting for GSE147520 and GSE207206.
- AnnData 2024 citation details, depending on target journal requirements.
- Qi et al. 2022 citation details.
- Whether the target journal prefers dataset citations in Methods/Data availability rather than numbered references.
"""

notes_out.write_text(notes)

print("Done.")
print(f"Wrote reference-checked manuscript: {out_md}")
print(f"Wrote checked reference table: {out_refs}")
print(f"Wrote verification notes: {notes_out}")
print()
print("Reference status counts:")
print(df["status"].value_counts().to_string())
