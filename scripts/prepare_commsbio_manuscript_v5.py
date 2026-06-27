#!/usr/bin/env python3

from pathlib import Path
import re

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"
DOC_DIR.mkdir(parents=True, exist_ok=True)

in_md = DOC_DIR / "manuscript_draft_v4_health_checked.md"
out_md = DOC_DIR / "manuscript_draft_v5_commsbio_ready.md"
checklist_out = DOC_DIR / "communications_biology_compliance_checklist_v1.md"
change_log = DOC_DIR / "manuscript_v5_commsbio_change_log.md"

if not in_md.exists():
    raise FileNotFoundError(in_md)

text = in_md.read_text()

# -------------------------
# Replace structured abstract with unstructured Communications Biology style
# -------------------------
new_abstract = """## Abstract

The human thymus supports T cell development through specialised epithelial, stromal and haematopoietic antigen presenting compartments, but the spatial organisation of stromal immune antigen presenting niches in early life remains incompletely resolved. Here we integrate public human thymus single cell and spatial transcriptomic datasets to reconstruct thymic epithelial, stromal and immune antigen presenting programmes. The Bautista human thymic stromal dataset was used as a reference for epithelial and stromal architecture, while the Heimli paediatric thymus CITE seq and Visium dataset was used to resolve immune antigen presenting, thymic B cell, epithelial, stromal and spatial niche programmes. Bautista analysis resolved broad TEC, fibroblast, vascular, vascular MHCII associated and perivascular compartments. Heimli analysis identified paediatric B cell APC, macrophage APC, activated DC like, pDC like, cTEC like, mTEC like, thymocyte, endothelial, fibroblast and perivascular states. Spatial projection across eight Heimli Visium sections identified recurrent MHCII APC rich regions that co patterned most strongly with B cell APC signatures, followed by mTEC like, macrophage APC and pDC like programmes, whereas cTEC like and T cell thymocyte signatures showed relative spatial separation. These findings support a multi component model in which human thymic antigen presenting niches are spatially patterned across epithelial and haematopoietic APC associated programmes, providing a reproducible framework for studying thymic stromal immune organisation."""

text = re.sub(
    r"## Abstract\n\n.*?\n\n## Introduction",
    new_abstract + "\n\n## Introduction",
    text,
    flags=re.S
)

# -------------------------
# Shorten Results headings to fit Communications Biology checklist
# -------------------------
heading_replacements = {
    "### Complementary thymus single-cell datasets support a two-arm strategy for reconstructing stromal and immune antigen-presenting niches":
        "### Two datasets define niche analysis arms",
    "### Heimli CITE-seq resolves B-cell, myeloid and plasmacytoid antigen-presenting compartments within paediatric thymus":
        "### Heimli resolves paediatric APC states",
    "### Spatial Visium projection identifies recurrent MHC-II-rich thymic regions associated with B-cell APC, mTEC-like and myeloid APC programmes":
        "### Visium maps recurrent MHCII rich regions",
    "### Integrated single-cell and spatial evidence supports a multi-component model of thymic antigen-presenting niches":
        "### Integrated evidence supports niche organisation",
}

for old, new in heading_replacements.items():
    text = text.replace(old, new)

# -------------------------
# Shorten Methods headings and remove punctuation/hyphens where possible
# -------------------------
methods_heading_replacements = {
    "### Single-cell preprocessing and annotation": "### Single cell preprocessing and annotation",
    "### Spatial co-pattern analysis": "### Spatial co pattern analysis",
}

for old, new in methods_heading_replacements.items():
    text = text.replace(old, new)

# -------------------------
# Move availability statements into End Notes style after References if needed
# Communications Biology order: Methods, References, End Notes, Figure legends.
# We keep Data/code in End Notes here for clarity; checklist accepts at end of Methods/main text.
# -------------------------
data_code_block = """## End Notes

### Data availability

All datasets analysed in this study are publicly available through the Gene Expression Omnibus. The Bautista human thymic stromal single cell dataset was obtained from GSE147520. The Heimli paediatric thymus CITE seq and Visium dataset was obtained from GSE207206. Processed analytical tables, figure source tables, annotation files and reproducibility outputs generated in this study are organised in the project directory and will be deposited in a public repository before publication.

### Code availability

Custom scripts used for data retrieval checks, quality control, single cell processing, RNA and ADT marker separation, spatial signature scoring, spatial co pattern analysis, figure generation and manuscript package generation are included in the project reproducibility folder and will be deposited in a public repository before publication.

### Author contributions

R.S.M. conceived the study, designed the computational workflow, performed data analysis, interpreted the results and drafted the manuscript. Additional author contributions will be added according to project participation.

### Acknowledgements

The author acknowledges the investigators who generated and made publicly available the human thymus single cell and spatial transcriptomic datasets re analysed in this study.

### Competing interests

The author declares no competing interests.
"""

# Remove old availability / end matter sections before References.
patterns_to_remove = [
    r"\n## Data and code availability\n\n.*?(?=\n## Author contributions)",
    r"\n## Author contributions\n\n.*?(?=\n## Acknowledgements)",
    r"\n## Acknowledgements\n\n.*?(?=\n## Competing interests)",
    r"\n## Competing interests\n\n.*?(?=\n## References)",
]
for pat in patterns_to_remove:
    text = re.sub(pat, "", text, flags=re.S)

# Insert End Notes after References block and before Figure legends, if Figure legends currently before References then reorder.
# Current structure likely has Figure legends before old endmatter and References; reconstruct from section blocks.

def get_section(src, heading):
    m = re.search(rf"^## {re.escape(heading)}\n", src, flags=re.M)
    if not m:
        return ""
    start = m.start()
    next_m = re.search(r"^## ", src[m.end():], flags=re.M)
    end = m.end() + next_m.start() if next_m else len(src)
    return src[start:end].strip()

title_part = text.split("\n## Abstract\n")[0].strip()
abstract = get_section(text, "Abstract")
intro = get_section(text, "Introduction")
results = get_section(text, "Results")
discussion = get_section(text, "Discussion")
methods = get_section(text, "Methods")
figlegs = get_section(text, "Figure legends")
refs = get_section(text, "References")

ordered = [
    title_part,
    abstract,
    intro,
    results,
    discussion,
    methods,
    refs,
    data_code_block.strip(),
    figlegs,
]

text_v5 = "\n\n".join([x for x in ordered if x.strip()]) + "\n"

# -------------------------
# Write output manuscript
# -------------------------
out_md.write_text(text_v5)

# -------------------------
# Compliance checks
# -------------------------
def section_word_count(src, heading):
    sec = get_section(src, heading)
    sec = re.sub(r"^## .+\n", "", sec)
    return len(re.findall(r"\b[\w'-]+\b", sec))

intro_wc = section_word_count(text_v5, "Introduction")
results_wc = section_word_count(text_v5, "Results")
discussion_wc = section_word_count(text_v5, "Discussion")
ird_wc = intro_wc + results_wc + discussion_wc
abstract_wc = section_word_count(text_v5, "Abstract")

result_headings = re.findall(r"^### (.+)$", get_section(text_v5, "Results"), flags=re.M)
methods_headings = re.findall(r"^### (.+)$", get_section(text_v5, "Methods"), flags=re.M)

def heading_ok(h):
    return len(h) < 60 and not any(ch in h for ch in [":", ";", ",", "?", "!", "(", ")", "/", "–", "-", "—"])

check_lines = []
check_lines.append("# Communications Biology compliance checklist v1\n")
check_lines.append(f"- Manuscript: `{out_md}`")
check_lines.append(f"- Abstract word count: {abstract_wc}")
check_lines.append(f"- Introduction word count: {intro_wc}")
check_lines.append(f"- Results word count: {results_wc}")
check_lines.append(f"- Discussion word count: {discussion_wc}")
check_lines.append(f"- Introduction + Results + Discussion word count: {ird_wc}")
check_lines.append("")
check_lines.append("## Checklist")
check_lines.append(f"- Main text order adapted to Title, Abstract, Introduction, Results, Discussion, Methods, References, End Notes, Figure legends: YES")
check_lines.append(f"- Introduction below 1,000 words: {'YES' if intro_wc < 1000 else 'NO'}")
check_lines.append(f"- Introduction + Results + Discussion below 5,000 words: {'YES' if ird_wc <= 5000 else 'NO'}")
check_lines.append(f"- Data availability included: {'YES' if '### Data availability' in text_v5 else 'NO'}")
check_lines.append(f"- Code availability included: {'YES' if '### Code availability' in text_v5 else 'NO'}")
check_lines.append(f"- Author contributions included: {'YES' if '### Author contributions' in text_v5 else 'NO'}")
check_lines.append(f"- Competing interests included: {'YES' if '### Competing interests' in text_v5 else 'NO'}")
check_lines.append("")
check_lines.append("## Results headings")
for h in result_headings:
    check_lines.append(f"- `{h}` | characters: {len(h)} | compliant: {'YES' if heading_ok(h) else 'NO'}")
check_lines.append("")
check_lines.append("## Methods headings")
for h in methods_headings:
    check_lines.append(f"- `{h}` | characters: {len(h)} | compliant: {'YES' if heading_ok(h) else 'NO'}")
checklist_out.write_text("\n".join(check_lines) + "\n")

# -------------------------
# Change log
# -------------------------
change_lines = []
change_lines.append("# Manuscript v5 Communications Biology change log\n")
change_lines.append(f"- Input: `{in_md}`")
change_lines.append(f"- Output: `{out_md}`")
change_lines.append("")
change_lines.append("## Changes")
change_lines.append("- Converted structured abstract into an unstructured abstract.")
change_lines.append("- Shortened Results headings to fewer than 60 characters without punctuation.")
change_lines.append("- Shortened Methods headings by removing hyphenated punctuation where possible.")
change_lines.append("- Reordered manuscript toward Communications Biology style.")
change_lines.append("- Added End Notes with Data availability, Code availability, Author contributions, Acknowledgements and Competing interests.")
change_lines.append("- Preserved cautious claim language around co patterning, deconvolution and direct cell contact.")
change_log.write_text("\n".join(change_lines) + "\n")

print("Done.")
print(f"Wrote Communications Biology manuscript: {out_md}")
print(f"Wrote compliance checklist: {checklist_out}")
print(f"Wrote change log: {change_log}")
print()
print(f"Abstract words: {abstract_wc}")
print(f"Introduction words: {intro_wc}")
print(f"Results words: {results_wc}")
print(f"Discussion words: {discussion_wc}")
print(f"Intro + Results + Discussion words: {ird_wc}")
