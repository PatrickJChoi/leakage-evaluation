import re
import os

preprint_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Leakage_Evaluation_Preprint.md")

with open(preprint_path, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Patterns to check
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
    "CLASSIC I": r"CLASSIC\s+I"
}

# Gene list to check italicization
genes_to_check = ["TNK1", "IGSF8", "KCNQ1", "CYP4F12", "SSTR1", "SIGMAR1", "DANT2", "NOXA1", 
                  "TBC1D10B", "TNRC6B", "FOXP1", "SELL", "TMEM63A", "GOLT1A"]

print("=== STARTING AUDIT ===")
issues_found = 0

for line_idx, line in enumerate(lines, 1):
    # Check invalid patterns
    for name, pattern in invalid_patterns.items():
        if re.search(pattern, line, re.IGNORECASE):
            print(f"Issue found: '{name}' on Line {line_idx}: {line.strip()[:100]}...")
            issues_found += 1
            
    # Check gene italicization
    for gene in genes_to_check:
        # Find occurrences of the gene name
        for match in re.finditer(r'\b' + gene + r'\b', line):
            # Check if it is enclosed in '*' (italicized)
            start, end = match.start(), match.end()
            # Look back and look ahead to see if they are asterisks
            is_italic = False
            if start > 0 and end < len(line):
                # Simple check for *Gene* or _Gene_
                before = line[max(0, start-2):start]
                after = line[end:min(len(line), end+2)]
                if '*' in before and '*' in after or '_' in before and '_' in after:
                    is_italic = True
            if not is_italic:
                # Exclude figure legends or table titles if they are part of a longer formatting, 
                # but standard practice is to italicize everywhere in text
                print(f"Issue found: Unitalicized gene '{gene}' on Line {line_idx} at char {start}: {line.strip()[max(0, start-20):end+20]}...")
                issues_found += 1

print(f"=== AUDIT COMPLETE: {issues_found} ISSUES FOUND ===")

# Calculate word count per section
print("\n=== SECTION WORD COUNTS ===")
current_section = "Header"
section_words = 0
for line in lines:
    if line.startswith("# "):
        # Print previous section
        if section_words > 0:
            print(f"{current_section}: {section_words} words")
        current_section = line.strip("#").strip()
        section_words = 0
    elif line.startswith("## "):
        if section_words > 0:
            print(f"{current_section}: {section_words} words")
        current_section = line.strip("#").strip()
        section_words = 0
    else:
        section_words += len(line.split())
if section_words > 0:
    print(f"{current_section}: {section_words} words")
