import urllib.request
import json
import gzip
import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from collections import Counter
import time

# ============================================================
# 1. Map gene symbols to Entrez Gene IDs via NCBI E-utilities
# ============================================================
print("=" * 60)
print("MAPPING GENE SYMBOLS TO ENTREZ IDs (NCBI E-utilities)")
print("=" * 60)

protein_coding_genes = [
    ("TBC1D10B", 18), ("TNRC6B", 17), ("FOXP1", 10),
    ("DHCR24", 3), ("G6PD", 2), ("SELL", 2),
    ("ADA", 1), ("RIC3", 1), ("EZH2", 1),
    ("UBE2M", 1), ("GTF2F1", 1), ("DCLRE1B", 1), ("BLMH", 1),
]

symbol_to_entrez = {}
for symbol, folds in protein_coding_genes:
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
           f"db=gene&term={symbol}[Gene+Name]+AND+Homo+sapiens[Organism]"
           f"&retmode=json&retmax=5")
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if id_list:
            entrez_id = int(id_list[0])
            symbol_to_entrez[symbol] = entrez_id
            print(f"  {symbol:<15} -> Entrez {entrez_id}")
        else:
            print(f"  {symbol:<15} -> NOT FOUND")
    except Exception as e:
        print(f"  {symbol:<15} -> ERROR: {e}")
    time.sleep(0.4)  # rate limit

print(f"\nMapped {len(symbol_to_entrez)}/{len(protein_coding_genes)} genes")

# ============================================================
# 2. Load expression matrix and metadata
# ============================================================
print("\n" + "=" * 60)
print("LOADING DATA")
print("=" * 60)

expr = pd.read_csv("GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz",
                    index_col=0, compression='gzip')
print(f"Expression matrix: {expr.shape}")

# Metadata
metadata = {}
with gzip.open("GSE91061_series_matrix.txt.gz", 'rt', encoding='utf-8', errors='ignore') as f:
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

# ============================================================
# 3. Filter: pre-treatment + binary response
# ============================================================
print("\n" + "=" * 60)
print("FILTERING: PRE-TREATMENT + BINARY RESPONSE")
print("=" * 60)

pre_meta = meta_df[meta_df['visit (pre or on treatment)'] == 'Pre'].copy()
print(f"Pre-treatment samples: {len(pre_meta)}")
print(f"Response breakdown: {dict(Counter(pre_meta['response']))}")

# Binarize: PRCR = Responder, PD = Non-responder, exclude SD/UNK
binary_meta = pre_meta[pre_meta['response'].isin(['PD', 'PRCR'])].copy()
binary_meta['binary_response'] = binary_meta['response'].map({
    'PRCR': 'Responder', 'PD': 'Non-responder'
})

print(f"\nAfter excluding SD and UNK:")
resp_counts = Counter(binary_meta['binary_response'])
print(f"  Total:           {len(binary_meta)}")
print(f"  Responders:      {resp_counts['Responder']}")
print(f"  Non-responders:  {resp_counts['Non-responder']}")

# Get matching expression columns
sample_ids = binary_meta['title'].tolist()
available = [s for s in sample_ids if s in expr.columns]
print(f"  In expression:   {len(available)}")

binary_meta = binary_meta[binary_meta['title'].isin(available)].set_index('title')

# ============================================================
# 4. Check which Entrez IDs are in the expression matrix
# ============================================================
print("\n" + "=" * 60)
print("GENE PRESENCE CHECK IN EXPRESSION MATRIX")
print("=" * 60)

expr_entrez_ids = set(expr.index.astype(int) if expr.index.dtype != object else 
                      [int(x) for x in expr.index])
print(f"Total genes in expression matrix: {len(expr_entrez_ids)}")

found_genes = []
missing_genes = []
for symbol, folds in protein_coding_genes:
    if symbol in symbol_to_entrez:
        eid = symbol_to_entrez[symbol]
        if eid in expr_entrez_ids:
            found_genes.append((symbol, eid, folds))
            print(f"  FOUND:   {symbol:<15} (Entrez {eid}, {folds}/18 folds)")
        else:
            missing_genes.append((symbol, eid, folds))
            print(f"  MISSING: {symbol:<15} (Entrez {eid} not in matrix)")
    else:
        print(f"  UNMAPPED: {symbol:<15}")

print(f"\nGenes found in GSE91061: {len(found_genes)}/{len(protein_coding_genes)}")

# ============================================================
# 5. Mann-Whitney U tests
# ============================================================
print("\n" + "=" * 60)
print("MANN-WHITNEY U TESTS: RESPONDER vs NON-RESPONDER")
print("=" * 60)

resp_samples = binary_meta.index[binary_meta['binary_response'] == 'Responder'].tolist()
nonresp_samples = binary_meta.index[binary_meta['binary_response'] == 'Non-responder'].tolist()
print(f"Responder samples:     {len(resp_samples)}")
print(f"Non-responder samples: {len(nonresp_samples)}")
print()

header = f"{'Gene':<15} {'Entrez':>8} {'Folds':>6} {'Mean_R':>10} {'Mean_NR':>10} {'Diff':>8} {'Direction':>10} {'U-stat':>8} {'p-value':>10} {'Sig':>4}"
print(header)
print("-" * len(header))

results = []
sig_count = 0

for symbol, eid, folds in found_genes:
    r_vals = expr.loc[eid, resp_samples].values.astype(float)
    nr_vals = expr.loc[eid, nonresp_samples].values.astype(float)
    
    mean_r = np.mean(r_vals)
    mean_nr = np.mean(nr_vals)
    diff = mean_r - mean_nr
    direction = "R > NR" if diff > 0 else "NR > R"
    
    stat, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    
    sig = "*" if pval < 0.05 else ""
    if pval < 0.05:
        sig_count += 1
    
    print(f"{symbol:<15} {eid:>8} {folds:>5}  {mean_r:>10.3f} {mean_nr:>10.3f} {diff:>+8.3f} {direction:>10} {stat:>8.0f} {pval:>10.4f} {sig:>4}")
    
    results.append({
        'gene': symbol, 'entrez_id': eid, 'folds_selected': folds,
        'mean_R': mean_r, 'mean_NR': mean_nr, 'diff_R_minus_NR': diff,
        'direction': direction, 'U_stat': stat, 'p_value': pval
    })

# ============================================================
# 6. Summary
# ============================================================
print(f"\n{'=' * 60}")
print(f"SUMMARY")
print(f"{'=' * 60}")
print(f"Protein-coding genes in top 20:            {len(protein_coding_genes)}")
print(f"Successfully mapped to Entrez IDs:         {len(symbol_to_entrez)}")
print(f"Found in GSE91061 expression matrix:       {len(found_genes)}")
print(f"Genes with p < 0.05 (Mann-Whitney U):      {sig_count}")
print(f"Genes with p >= 0.05:                      {len(found_genes) - sig_count}")

# Direction consistency check
print(f"\nDirection of effect (in rld space):")
for r in results:
    consistency = ""
    # In our scRNA-seq model, higher expression in R was the training signal
    print(f"  {r['gene']:<15}  {r['direction']:<10}  p={r['p_value']:.4f}  "
          f"diff={r['diff_R_minus_NR']:+.3f}")

# Save
results_df = pd.DataFrame(results)
results_df.to_csv("gse91061_validation_results.csv", index=False)
print(f"\nResults saved to gse91061_validation_results.csv")
