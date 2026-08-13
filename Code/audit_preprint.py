import re
import os
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Manuscripts checked by this audit. Both the Markdown preprint and the
# LaTeX/BMC build must be checked -- numbers have drifted between the two
# before, and a check against only one file can silently miss the other.
# filetype is declared explicitly here ("md" / "tex") rather than inferred
# from content, so downstream checks select their rule by file type, not
# by guessing from what a given line looks like.
TARGET_FILES = [
    ("Leakage_Evaluation_Preprint.md", os.path.join(BASE_DIR, "Leakage_Evaluation_Preprint.md"), "md"),
    ("latex/leakage_evaluation_bmc.tex", os.path.join(BASE_DIR, "latex", "leakage_evaluation_bmc.tex"), "tex"),
    ("README.md", os.path.join(BASE_DIR, "README.md"), "md"),
]

missing = [label for label, path, _ in TARGET_FILES if not os.path.isfile(path)]
if missing:
    for label in missing:
        print(f"ERROR: required manuscript file not found: {label}")
    sys.exit(2)

file_lines = {}
file_types = {}
for label, path, filetype in TARGET_FILES:
    with open(path, "r", encoding="utf-8") as f:
        file_lines[label] = f.readlines()
    file_types[label] = filetype

# Patterns to check
# --------------------------------------------------------------------------
# RETIRED 2026-07-28: the bulk-validation patterns below (0/13, 0 of 13,
# 10/13, 10 of 13, 10 out of 13, 76.9, 0.092, 13 protein-coding) were added
# when the melanoma external bulk-validation gene set moved from an
# undocumented protein-coding filter (13 genes) to a stated rule -- a gene
# is tested if it was selected in >= 2 of the 18 melanoma LOO-CV folds AND
# is testable (resolves to an Entrez ID present in the GSE91061 matrix) --
# which yields 6 genes. The superseding values are 0/6, 5/6, 83.3%, and
# 0.219. Any reappearance of the values below is a regression back to the
# retired numbers, not a legitimate new use of those strings.
# --------------------------------------------------------------------------
invalid_patterns = {
    "75.0%": r"75\.0%",
    "p=0.0240": r"p\s*=\s*0\.0240",
    "FAM83E": r"\bFAM83E\b",
    "ARHGEF38": r"\bARHGEF38\b",
    "CSKMT": r"\bCSKMT\b",
    "MLXIPL": r"\bMLXIPL\b",
    "[cite] placeholder": r"\[cite\]",
    "TAURUS Consortium 2024": r"TAURUS\s+Consortium\s+2024",
    "TAURUS 1780-1794": r"1780\s*-\s*1794",
    "West et al. 2017": r"West\s+et\s+al",
    "Oncostatin M": r"Oncostatin\s+M",
    "Hanauer 2006": r"Hanauer",
    "CLASSIC I": r"CLASSIC\s+I",
    "0/13": r"0/13",
    "0 of 13": r"0\s+of\s+13",
    "10/13": r"10/13",
    "10 of 13": r"10\s+of\s+13",
    "10 out of 13": r"10\s+out\s+of\s+13",
    "76.9": r"76\.9",
    "0.092": r"0\.092",
    "13 protein-coding": r"13\s+protein-coding",
}

# Gene list to check italicization
# NOTE (2026-07-28): does not include DHCR24 or G6PD. Flagged for review,
# not added automatically -- confirm intended usage before extending this
# list.
genes_to_check = ["TNK1", "IGSF8", "KCNQ1", "CYP4F12", "SSTR1", "SIGMAR1", "DANT2", "NOXA1",
                  "TBC1D10B", "TNRC6B", "FOXP1", "SELL", "TMEM63A", "GOLT1A"]

TEX_ITALIC_RE = re.compile(r'\\(?:textit|emph)\{([^}]*)\}')


def tex_italic_spans(line):
    """Character spans (start, end) of text inside \\textit{...} / \\emph{...} on this line."""
    return [(m.start(1), m.end(1)) for m in TEX_ITALIC_RE.finditer(line)]


def is_italicized_tex(spans, start, end):
    return any(start >= s and end <= e for s, e in spans)


def is_italicized_md(line, start, end):
    # Simple check for *Gene* or _Gene_
    before = line[max(0, start - 2):start]
    after = line[end:min(len(line), end + 2)]
    return ('*' in before and '*' in after) or ('_' in before and '_' in after)


# References-section exclusion for the gene-italicization check ONLY.
# Reference-list titles reproduce a cited paper's own title text verbatim
# (e.g. "...kinase 1 (TNK1) facilitates..."), which is never italicized by
# citation convention regardless of whether TNK1 is also a tracked gene --
# so this is not a real formatting defect and should not be flagged. The
# forbidden-pattern check above is NOT scoped by this exclusion; it still
# runs over every line of the file.
# md: every line from (and including) the "# References" heading up to
# (but not including) the next top-level "# " heading -- i.e. just the
# References section, not everything after it.
# tex: lines inside a \begin{thebibliography}...\end{thebibliography}
# environment, if the source file contains one. As of 2026-07-28 the
# committed .tex does NOT contain such an environment -- references are
# compiled separately by bibtex8 from latex/references.bib into the .bbl,
# which this script does not read -- so this exclusion is currently a
# no-op for the .tex target and is included only so the check is correct
# if a literal bibliography is ever pasted into the .tex source.
gene_check_excluded = {}
for label, lines in file_lines.items():
    filetype = file_types[label]
    excluded = set()
    if filetype == "md":
        in_references = False
        for line_idx, line in enumerate(lines, 1):
            stripped = line.strip()
            if re.match(r'^#\s*References\s*$', stripped):
                in_references = True
                excluded.add(line_idx)
                continue
            if in_references:
                if re.match(r'^#\s', stripped):
                    in_references = False
                else:
                    excluded.add(line_idx)
    else:  # tex
        in_bibliography = False
        for line_idx, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith(r'\begin{thebibliography}'):
                in_bibliography = True
            if in_bibliography:
                excluded.add(line_idx)
            if stripped.startswith(r'\end{thebibliography}'):
                in_bibliography = False
    gene_check_excluded[label] = excluded

print("=== STARTING AUDIT ===")
issues_found = 0

for label, lines in file_lines.items():
    filetype = file_types[label]
    for line_idx, line in enumerate(lines, 1):
        # Check invalid patterns -- runs over the whole file, unscoped by
        # the References exclusion above.
        for name, pattern in invalid_patterns.items():
            if re.search(pattern, line, re.IGNORECASE):
                print(f"Issue found [{label}]: '{name}' on Line {line_idx}: {line.strip()[:100]}...")
                issues_found += 1

        if line_idx in gene_check_excluded[label]:
            continue

        # Check gene italicization -- rule selected by declared file type
        if filetype == "tex":
            spans = tex_italic_spans(line)
        else:
            spans = None
        for gene in genes_to_check:
            # Find occurrences of the gene name
            for match in re.finditer(r'\b' + gene + r'\b', line):
                start, end = match.start(), match.end()
                if filetype == "tex":
                    is_italic = is_italicized_tex(spans, start, end)
                else:
                    is_italic = is_italicized_md(line, start, end)
                if not is_italic:
                    # Exclude figure legends or table titles if they are part of a longer formatting,
                    # but standard practice is to italicize everywhere in text
                    print(f"Issue found [{label}]: Unitalicized gene '{gene}' on Line {line_idx} at char {start}: {line.strip()[max(0, start-20):end+20]}...")
                    issues_found += 1

print(f"=== AUDIT COMPLETE: {issues_found} ISSUES FOUND ===")

# Calculate word count per section -- section-boundary rule selected by
# declared file type, not guessed from line content.
TEX_SECTION_RE = re.compile(r'^\\(?:sub)?section\*\{(.*)\}\s*$')

for label, lines in file_lines.items():
    filetype = file_types[label]
    print(f"\n=== SECTION WORD COUNTS [{label}] ===")
    current_section = "Header"
    section_words = 0
    for line in lines:
        if filetype == "md":
            if line.startswith("# ") or line.startswith("## "):
                if section_words > 0:
                    print(f"{current_section}: {section_words} words")
                current_section = line.strip("#").strip()
                section_words = 0
            else:
                section_words += len(line.split())
        else:  # tex
            m = TEX_SECTION_RE.match(line.strip())
            if m:
                if section_words > 0:
                    print(f"{current_section}: {section_words} words")
                current_section = m.group(1).strip()
                section_words = 0
            else:
                section_words += len(line.split())
    if section_words > 0:
        print(f"{current_section}: {section_words} words")

sys.exit(1 if issues_found > 0 else 0)
