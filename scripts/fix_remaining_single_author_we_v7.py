#!/usr/bin/env python3

from pathlib import Path
import re
import shutil

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"

in_md = DOC_DIR / "manuscript_draft_v6_single_author_commsbio.md"
out_md = DOC_DIR / "manuscript_draft_v7_single_author_clean_commsbio.md"
report = DOC_DIR / "manuscript_v7_single_author_clean_report.md"

if not in_md.exists():
    raise FileNotFoundError(in_md)

text = in_md.read_text()

replacements = {
    "Importantly, we retained mixed or low-confidence clusters with conservative labels and avoided using these regions as the basis for strong lineage-specific claims.":
    "Importantly, mixed or low-confidence clusters were retained with conservative labels and were not used as the basis for strong lineage-specific claims.",
}

changes = []
for old, new in replacements.items():
    count = text.count(old)
    if count:
        text = text.replace(old, new)
        changes.append((old, new, count))

out_md.write_text(text)

remaining = []
for pat in [r"\bwe\b", r"\bWe\b", r"\bour\b", r"\bOur\b", r"\bus\b", r"\bUs\b"]:
    for m in re.finditer(pat, text):
        snippet = text[max(0, m.start()-80):m.end()+80].replace("\n", " ")
        remaining.append((pat, snippet))

lines = []
lines.append("# Manuscript v7 single-author clean report\n")
lines.append(f"- Input: `{in_md}`")
lines.append(f"- Output: `{out_md}`")
lines.append("")
lines.append("## Changes made")
if changes:
    for old, new, count in changes:
        lines.append(f"- Replaced sentence: {count} occurrence(s).")
else:
    lines.append("- No target sentence was found.")
lines.append("")
lines.append("## Remaining first-person plural wording")
if remaining:
    for pat, snippet in remaining:
        lines.append(f"- `{pat}`: {snippet}")
else:
    lines.append("- No remaining `we/our/us` wording detected.")

report.write_text("\n".join(lines) + "\n")

print("Done.")
print(f"Wrote clean manuscript: {out_md}")
print(f"Wrote report: {report}")
print(f"Sentence replacements: {sum(c for _, _, c in changes)}")
print(f"Remaining plural first-person hits: {len(remaining)}")
