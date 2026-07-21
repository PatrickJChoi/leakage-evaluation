"""
Reproduce the GSE91061 directional-concordance count and binomial test
end-to-end from raw inputs, without hardcoding n_agree=10.

Reuses the IDENTICAL raw-file parsing and direction logic as
validate_gse91061.py (metadata parsing: lines ~58-119 there; direction rule:
lines 172-175 there). Not imported directly because validate_gse91061.py
performs live NCBI E-utilities calls at module level on import; instead the
same lines are reproduced verbatim here for the metadata/expression parsing
and the direction computation, so direction cannot drift from the values
computed in validate_gse91061.py. Entrez gene IDs are taken from validate_gse91061_fpkm.py:10-13
(the already-resolved IDs for this same frozen 13-gene signature), avoiding
a repeat live NCBI dependency.

Honest recompute; targets no value.
"""
import gzip, json
import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import mannwhitneyu, binomtest

RAW_EXPR_FILE = "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"
RAW_META_FILE = "GSE91061_series_matrix.txt.gz"

# Frozen 13-gene discovery signature (symbol, entrez_id, folds_selected),
# identical set as validate_gse91061.py:18-23; entrez IDs as already resolved
# in validate_gse91061_fpkm.py:10-13.
genes = [
    ("TBC1D10B", 26000, 18), ("TNRC6B", 23112, 17), ("FOXP1", 27086, 10),
    ("DHCR24", 1718, 3), ("G6PD", 2539, 2), ("SELL", 6402, 2),
    ("ADA", 100, 1), ("RIC3", 79608, 1), ("EZH2", 2146, 1),
    ("UBE2M", 9040, 1), ("GTF2F1", 2962, 1), ("DCLRE1B", 64858, 1),
    ("BLMH", 642, 1),
]

# --- 1. Load expression matrix (identical to validate_gse91061.py:54-56) ---
expr = pd.read_csv(RAW_EXPR_FILE, index_col=0, compression='gzip')
print(f"Expression matrix: {expr.shape}")

# --- Metadata parsing (identical to validate_gse91061.py:59-89) ---
metadata = {}
with gzip.open(RAW_META_FILE, 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        ls = line.strip()
        if ls.startswith('!'):
            parts = ls.split('\t')
            key = parts[0]
            vals = [v.strip('"') for v in parts[1:]]
            if key not in metadata:
                metadata[key] = []
            metadata[key].append(vals)
        elif ls.startswith('"ID_REF"'):
            break

geo_ids = metadata["!Sample_geo_accession"][0]
titles = metadata["!Sample_title"][0]
chars = metadata.get("!Sample_characteristics_ch1", [])

sample_meta = []
for i in range(len(geo_ids)):
    rec = {"geo": geo_ids[i], "title": titles[i]}
    for char_line in chars:
        if i < len(char_line):
            val = char_line[i]
            if ':' in val:
                k = val.split(':')[0].strip().lower()
                v = val.split(':', 1)[1].strip()
                rec[k] = v
    sample_meta.append(rec)

meta_df = pd.DataFrame(sample_meta)

# --- Filter: pre-treatment + binary response (identical to lines 98-119) ---
pre_meta = meta_df[meta_df['visit (pre or on treatment)'] == 'Pre'].copy()
print(f"Pre-treatment samples: {len(pre_meta)}")
print(f"Response breakdown: {dict(Counter(pre_meta['response']))}")

binary_meta = pre_meta[pre_meta['response'].isin(['PD', 'PRCR'])].copy()
binary_meta['binary_response'] = binary_meta['response'].map({
    'PRCR': 'Responder', 'PD': 'Non-responder'
})

resp_counts = Counter(binary_meta['binary_response'])
print(f"After excluding SD and UNK: total={len(binary_meta)} "
      f"Responders={resp_counts['Responder']} Non-responders={resp_counts['Non-responder']}")

sample_ids = binary_meta['title'].tolist()
available = [s for s in sample_ids if s in expr.columns]
print(f"In expression: {len(available)}")

binary_meta = binary_meta[binary_meta['title'].isin(available)].set_index('title')

resp_samples = binary_meta.index[binary_meta['binary_response'] == 'Responder'].tolist()
nonresp_samples = binary_meta.index[binary_meta['binary_response'] == 'Non-responder'].tolist()
print(f"Responder samples: {len(resp_samples)}  Non-responder samples: {len(nonresp_samples)}")

# --- 2/3. Per-gene direction (identical rule to lines 168-177) + tally ---
print(f"\n{'Gene':<12}{'Entrez':>8}{'mean_R':>10}{'mean_NR':>10}{'diff':>9}{'direction':>10}{'p_value':>10}")
rows = []
n_agree = 0  # count of "R > NR"
for symbol, eid, folds in genes:
    if eid not in expr.index:
        print(f"{symbol:<12} MISSING from expression matrix (Entrez {eid})")
        continue
    r_vals = expr.loc[eid, resp_samples].values.astype(float)
    nr_vals = expr.loc[eid, nonresp_samples].values.astype(float)

    mean_r = np.mean(r_vals)
    mean_nr = np.mean(nr_vals)
    diff = mean_r - mean_nr
    direction = "R > NR" if diff > 0 else "NR > R"
    if direction == "R > NR":
        n_agree += 1

    stat, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

    print(f"{symbol:<12}{eid:>8}{mean_r:>10.3f}{mean_nr:>10.3f}{diff:>+9.3f}{direction:>10}{pval:>10.4f}")
    rows.append({
        "gene": symbol, "entrez_id": eid, "folds_selected": folds,
        "mean_R": mean_r, "mean_NR": mean_nr, "diff_R_minus_NR": diff,
        "direction": direction, "p_value": pval,
    })

n_total = len(rows)
print(f"\nn_agree (R > NR) = {n_agree} / {n_total}")

# --- 4. Binomial tests, both alternatives, fed with the DERIVED n_agree ---
bt_greater = binomtest(n_agree, n_total, 0.5, alternative='greater')
bt_two_sided = binomtest(n_agree, n_total, 0.5, alternative='two-sided')
print(f"binomtest(n_agree={n_agree}, n={n_total}, p=0.5, alternative='greater')   = {bt_greater.pvalue}")
print(f"binomtest(n_agree={n_agree}, n={n_total}, p=0.5, alternative='two-sided') = {bt_two_sided.pvalue}")

# --- 5. Save everything ---
out = {
    "raw_inputs_used": [RAW_EXPR_FILE, RAW_META_FILE],
    "per_gene": rows,
    "n_agree": n_agree,
    "n_total": n_total,
    "binomial_p_one_sided_greater": bt_greater.pvalue,
    "binomial_p_two_sided": bt_two_sided.pvalue,
}
with open("reproduce_directional_binomial_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved reproduce_directional_binomial_results.json")
