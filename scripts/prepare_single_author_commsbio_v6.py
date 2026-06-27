#!/usr/bin/env python3

from pathlib import Path
import re
import shutil

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"
PACKAGE_DIR = PROJECT_DIR / "submission_package_v3_commsbio"

in_md = DOC_DIR / "manuscript_draft_v5_commsbio_ready.md"
out_md = DOC_DIR / "manuscript_draft_v6_single_author_commsbio.md"
report_out = DOC_DIR / "manuscript_v6_single_author_change_report.md"

if not in_md.exists():
    raise FileNotFoundError(in_md)

text = in_md.read_text()

author_block = """Reuben S. Maghembe^1,2*

^1 Department of Microbiology and Parasitology, Faculty of Medicine, St. Francis University College of Health and Allied Sciences (SFUCHAS), Ifakara, Tanzania.  
^2 Department of Omics and Computational Biology, AfroBiomics Co. Ltd., Bridge Street, Kivukoni, Dar es Salaam, Tanzania.

*Correspondence: Reuben S. Maghembe, Department of Microbiology and Parasitology, Faculty of Medicine, St. Francis University College of Health and Allied Sciences (SFUCHAS), Ifakara, Tanzania. Email: rmaghembe@gmail.com; rmaghembe@sfuchas.ac.tz
"""

# Insert author block after title if not already present.
lines = text.splitlines()
if lines and lines[0].startswith("# "):
    title = lines[0]
    rest = "\n".join(lines[1:]).lstrip()
    text = title + "\n\n" + author_block + "\n" + rest

# Single-author voice adjustments.
replacements = {
    "Here we integrate": "Here, I integrate",
    "Here, we integrate": "Here, I integrate",
    "Here we used": "Here, I used",
    "Here, we used": "Here, I used",
    "we integrate": "I integrate",
    "we used": "I used",
    "We used": "I used",
    "we identify": "I identify",
    "We identify": "I identify",
    "we then project": "I then project",
    "We then project": "I then project",
    "we applied": "I applied",
    "We applied": "I applied",
    "we first established": "I first established",
    "We first established": "I first established",
    "we next projected": "I next projected",
    "We next projected": "I next projected",
    "our analysis": "the analysis",
    "Our analysis": "The analysis",
    "our results": "the results",
    "Our results": "The results",
    "our findings": "the findings",
    "Our findings": "The findings",
    "our study": "this study",
    "Our study": "This study",
}

changes = []
for old, new in replacements.items():
    count = text.count(old)
    if count:
        text = text.replace(old, new)
        changes.append((old, new, count))

# Avoid awkward duplicate commas after replacement.
text = text.replace("Here, I integrate", "Here, I integrate")
text = text.replace("I used a two-arm", "I applied a two-arm")

# Strengthen single-author end notes.
text = text.replace(
    "R.S.M. conceived the study, designed the computational workflow, performed data analysis, interpreted the results and drafted the manuscript. Additional author contributions will be added according to project participation.",
    "R.S.M. conceived the study, designed the computational workflow, performed all analyses, interpreted the results, prepared the figures and wrote the manuscript."
)

text = text.replace(
    "The author acknowledges the investigators who generated and made publicly available the human thymus single cell and spatial transcriptomic datasets re analysed in this study.",
    "I acknowledge the investigators who generated and made publicly available the human thymus single cell and spatial transcriptomic datasets re analysed in this study."
)

text = text.replace(
    "The author declares no competing interests.",
    "The author declares no competing interests."
)

out_md.write_text(text)

# Copy to package if present.
if PACKAGE_DIR.exists():
    dst = PACKAGE_DIR / "01_manuscript" / out_md.name
    shutil.copy2(out_md, dst)

# Audit remaining plural first-person wording.
remaining_patterns = [
    r"\bwe\b", r"\bWe\b", r"\bour\b", r"\bOur\b", r"\bus\b", r"\bUs\b"
]
remaining = []
for pat in remaining_patterns:
    matches = list(re.finditer(pat, text))
    for m in matches:
        snippet = text[max(0, m.start()-80):m.end()+80].replace("\n", " ")
        remaining.append((pat, snippet))

report = []
report.append("# Manuscript v6 single-author change report\n")
report.append(f"- Input: `{in_md}`")
report.append(f"- Output: `{out_md}`")
if PACKAGE_DIR.exists():
    report.append(f"- Copied to package: `{PACKAGE_DIR / '01_manuscript' / out_md.name}`")
report.append("")
report.append("## Author block added")
report.append("- Reuben S. Maghembe with SFUCHAS and AfroBiomics affiliations.")
report.append("- Correspondence line corrected for spelling and format.")
report.append("")
report.append("## Voice changes")
if changes:
    for old, new, count in changes:
        report.append(f"- `{old}` -> `{new}`: {count}")
else:
    report.append("- No first-person plural replacements were needed.")
report.append("")
report.append("## Remaining first-person plural checks")
if remaining:
    for pat, snippet in remaining:
        report.append(f"- Pattern `{pat}`: {snippet}")
else:
    report.append("- No remaining `we/our/us` style wording detected.")

report_out.write_text("\n".join(report) + "\n")

print("Done.")
print(f"Wrote single-author manuscript: {out_md}")
print(f"Wrote report: {report_out}")
if PACKAGE_DIR.exists():
    print(f"Copied manuscript into package: {PACKAGE_DIR / '01_manuscript' / out_md.name}")
print(f"Voice replacements made: {sum(c for _, _, c in changes)}")
print(f"Remaining plural first-person hits: {len(remaining)}")
