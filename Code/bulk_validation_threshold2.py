"""
GSE91061 bulk validation of the melanoma consensus signature under a stated,
principled threshold: genes selected in >= 2 of the 18 melanoma LOO-CV folds
(replacing the prior undocumented top-20-then-manual-protein-coding cut).

The gene set is loaded programmatically from melanoma_consensus_counts.pkl,
a gene->fold-count mapping (no gene names are hardcoded here). Entrez Gene
IDs are loaded from the versioned mapping artifact Results/entrez_mapping.json
(produced in a prior read-only NCBI E-utilities session) rather than from a
hardcoded dict or a live network call — this script makes NO network calls
at all. melanoma_consensus_counts.pkl is tracked in git, whereas the
previously-used pickle (melanoma_verified_results.pkl) was not.

Explicit, stated exclusion rule: a gene is testable only if (1) it resolves
to a non-null Entrez ID in entrez_mapping.json, AND (2) that Entrez ID is
present in the GSE91061 expression matrix's index (checked directly against
the loaded matrix, not merely against the mapping file's cached flag). Every
excluded gene, its fold count, and which of the two conditions it failed are
reported to stdout and recorded in the output JSON.

Metadata parsing, cohort filtering, and direction logic are reproduced
identically to Code/reproduce_directional_binomial.py, for the same reason
that script gives for not importing validate_gse91061.py directly: avoiding
a live-network-call side effect on import (moot here since this script has
no network dependency at all, but the parsing/filtering/direction code is
kept byte-identical to avoid any drift from the values computed there).

Honest recompute; targets no value.
"""
import gzip, json, os, pickle
import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import mannwhitneyu, binomtest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # papers/leakage-evaluation
MELANOMA_COUNTS_PKL = os.path.join(REPO_ROOT, "melanoma_consensus_counts.pkl")
ENTREZ_MAPPING_JSON = os.path.join(REPO_ROOT, "Results", "entrez_mapping.json")
RAW_EXPR_FILE = os.path.join(REPO_ROOT, "Data", "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz")
RAW_META_FILE = os.path.join(REPO_ROOT, "Data", "GSE91061_series_matrix.txt.gz")
OUT_PATH = os.path.join(REPO_ROOT, "Results", "bulk_validation_threshold2.json")

FOLD_THRESHOLD = 2  # rule: genes selected in >= 2 of 18 melanoma LOO-CV folds

# --- 1. Load the >=2-fold gene set from melanoma_consensus_counts.pkl (not hardcoded) ---
with open(MELANOMA_COUNTS_PKL, "rb") as f:
    feature_freq = pickle.load(f)
gene_set = sorted(
    [(g, c) for g, c in feature_freq.items() if c >= FOLD_THRESHOLD],
    key=lambda gc: (-gc[1], gc[0]),
)
print(f"Genes selected in >= {FOLD_THRESHOLD} of 18 melanoma LOO-CV folds: {len(gene_set)}")
for gene, cnt in gene_set:
    print(f"  {gene}: {cnt}/18")

# --- 2. Load Entrez IDs from the versioned mapping artifact (no network calls) ---
with open(ENTREZ_MAPPING_JSON, "r", encoding="utf-8") as f:
    entrez_mapping = json.load(f)["genes"]

# --- 3. Load expression matrix and determine testability against it directly ---
expr = pd.read_csv(RAW_EXPR_FILE, index_col=0, compression='gzip')
print(f"\nExpression matrix: {expr.shape}")
expr_entrez_ids = set(
    expr.index.astype(int) if expr.index.dtype != object else [int(x) for x in expr.index]
)

# --- 4. Explicit, stated exclusion: testable iff (a) resolves to a non-null Entrez ID
#        in entrez_mapping.json AND (b) that ID is present in the GSE91061 matrix index.
testable = []
excluded = []
for gene, cnt in gene_set:
    entry = entrez_mapping.get(gene)
    if entry is None or entry.get("entrez_id") is None:
        excluded.append({
            "gene": gene, "folds_selected": cnt,
            "failed_condition": "no_entrez_id",
            "detail": "no non-null Entrez ID for this gene in Results/entrez_mapping.json",
        })
        continue
    eid = entry["entrez_id"]
    if eid not in expr_entrez_ids:
        excluded.append({
            "gene": gene, "folds_selected": cnt, "entrez_id": eid,
            "failed_condition": "entrez_id_not_in_matrix",
            "detail": f"Entrez ID {eid} resolved but not present in the GSE91061 expression matrix index",
        })
        continue
    testable.append((gene, eid, cnt))

print(f"\nTestable genes (resolved Entrez ID AND present in GSE91061 matrix): {len(testable)}")
for gene, eid, cnt in testable:
    print(f"  {gene}: entrez={eid}  folds={cnt}/18")

print(f"\nExcluded genes: {len(excluded)}")
for e in excluded:
    if e["failed_condition"] == "no_entrez_id":
        print(f"  {e['gene']}: folds={e['folds_selected']}/18  FAILED: no_entrez_id  ({e['detail']})")
    else:
        print(f"  {e['gene']}: folds={e['folds_selected']}/18  entrez={e['entrez_id']}  "
              f"FAILED: entrez_id_not_in_matrix  ({e['detail']})")

genes = testable

# --- Metadata parsing (identical to reproduce_directional_binomial.py) ---
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

# --- Filter: pre-treatment + binary response (identical to reproduce_directional_binomial.py) ---
pre_meta = meta_df[meta_df['visit (pre or on treatment)'] == 'Pre'].copy()
print(f"\nPre-treatment samples: {len(pre_meta)}")
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

# --- Per-gene direction (identical rule to reproduce_directional_binomial.py) + tally ---
print(f"\n{'Gene':<15}{'Entrez':>8}{'mean_R':>10}{'mean_NR':>10}{'diff':>9}{'direction':>10}{'p_value':>10}")
rows = []
n_agree = 0  # count of "R > NR"
for symbol, eid, folds in genes:
    if eid not in expr.index:
        print(f"{symbol:<15} MISSING from expression matrix (Entrez {eid})")
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

    print(f"{symbol:<15}{eid:>8}{mean_r:>10.3f}{mean_nr:>10.3f}{diff:>+9.3f}{direction:>10}{pval:>10.4f}")
    rows.append({
        "gene": symbol, "entrez_id": eid, "folds_selected": folds,
        "mean_R": mean_r, "mean_NR": mean_nr, "diff_R_minus_NR": diff,
        "direction": direction, "p_value": pval,
    })

n_total = len(rows)
print(f"\nn_agree (R > NR) = {n_agree} / {n_total}")

# --- Binomial tests, both alternatives ---
bt_greater = binomtest(n_agree, n_total, 0.5, alternative='greater')
bt_two_sided = binomtest(n_agree, n_total, 0.5, alternative='two-sided')
print(f"binomtest(n_agree={n_agree}, n={n_total}, p=0.5, alternative='greater')   = {bt_greater.pvalue}")
print(f"binomtest(n_agree={n_agree}, n={n_total}, p=0.5, alternative='two-sided') = {bt_two_sided.pvalue}")

# --- Save everything, including the exclusion record, to Results/ (tracked), not Data/ (gitignored) ---
out = {
    "threshold_rule": f">= {FOLD_THRESHOLD} of 18 melanoma LOO-CV folds",
    "testability_rule": "testable iff resolved to a non-null Entrez ID in Results/entrez_mapping.json AND that ID is present in the GSE91061 expression matrix index",
    "raw_inputs_used": [RAW_EXPR_FILE, RAW_META_FILE],
    "gene_source": MELANOMA_COUNTS_PKL,
    "entrez_mapping_source": ENTREZ_MAPPING_JSON,
    "excluded_genes": excluded,
    "per_gene": rows,
    "n_agree": n_agree,
    "n_total": n_total,
    "binomial_p_one_sided_greater": bt_greater.pvalue,
    "binomial_p_two_sided": bt_two_sided.pvalue,
}
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved {OUT_PATH}")
