"""
Re-derive GSE91061 bulk validation from the RAW FPKM file (not the saved CSV,
which was built from the rld/DESeq2 file). Entrez IDs taken from the existing
result table to avoid a live NCBI dependency. Honest recompute.
"""
import gzip, pandas as pd, numpy as np
from collections import Counter
from scipy.stats import mannwhitneyu

genes = [("TBC1D10B",26000,18),("TNRC6B",23112,17),("FOXP1",27086,10),
         ("DHCR24",1718,3),("G6PD",2539,2),("SELL",6402,2),("ADA",100,1),
         ("RIC3",79608,1),("EZH2",2146,1),("UBE2M",9040,1),("GTF2F1",2962,1),
         ("DCLRE1B",64858,1),("BLMH",642,1)]

expr = pd.read_csv("GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
                   index_col=0, compression='gzip')
expr.index = [int(x) for x in expr.index]
print(f"FPKM matrix: {expr.shape}")

# series-matrix metadata
metadata = {}
with gzip.open("GSE91061_series_matrix.txt.gz",'rt',encoding='utf-8',errors='ignore') as f:
    for line in f:
        ls = line.strip()
        if ls.startswith('!'):
            parts = ls.split('\t'); key = parts[0]
            metadata.setdefault(key, []).append([v.strip('"') for v in parts[1:]])
        elif ls.startswith('"ID_REF"'): break
titles = metadata["!Sample_title"][0]
chars = metadata.get("!Sample_characteristics_ch1", [])
recs = []
for i in range(len(titles)):
    rec = {"title": titles[i]}
    for cl in chars:
        if i < len(cl) and ':' in cl[i]:
            k = cl[i].split(':')[0].strip().lower(); v = cl[i].split(':',1)[1].strip()
            rec[k] = v
    recs.append(rec)
meta_df = pd.DataFrame(recs)

pre = meta_df[meta_df['visit (pre or on treatment)'] == 'Pre'].copy()
print(f"Pre-treatment samples: {len(pre)}  response breakdown: {dict(Counter(pre['response']))}")
binm = pre[pre['response'].isin(['PD','PRCR'])].copy()
binm['grp'] = binm['response'].map({'PRCR':'Responder','PD':'Non-responder'})
avail = [s for s in binm['title'] if s in expr.columns]
binm = binm[binm['title'].isin(avail)].set_index('title')
resp = binm.index[binm['grp']=='Responder'].tolist()
nonr = binm.index[binm['grp']=='Non-responder'].tolist()
print(f"Baseline used: total={len(binm)}  Responders(CR/PR)={len(resp)}  Non-responders(PD)={len(nonr)}")

print(f"\n{'Gene':<10}{'Entrez':>8}{'Folds':>6}{'meanR':>9}{'meanNR':>9}{'dir':>8}{'U':>7}{'p_fpkm':>9}")
rows = []; sig = 0; correct_dir = 0
for sym, eid, folds in genes:
    if eid not in expr.index:
        print(f"{sym:<10}{eid:>8}  NOT IN FPKM MATRIX"); continue
    rv = expr.loc[eid, resp].values.astype(float)
    nv = expr.loc[eid, nonr].values.astype(float)
    mr, mn = np.mean(rv), np.mean(nv); diff = mr-mn
    d = "R>NR" if diff>0 else "NR>R"
    if diff>0: correct_dir += 1
    U, p = mannwhitneyu(rv, nv, alternative='two-sided')
    if p<0.05: sig += 1
    print(f"{sym:<10}{eid:>8}{folds:>6}{mr:>9.3f}{mn:>9.3f}{d:>8}{U:>7.0f}{p:>9.4f}")
    rows.append({"gene":sym,"entrez":eid,"folds":folds,"meanR":mr,"meanNR":mn,
                 "diff":diff,"direction":d,"U":U,"p_value":p})
nfound = len(rows)
print(f"\nGenes found in FPKM: {nfound}/13 | p<0.05: {sig}/{nfound} | R>NR direction: {correct_dir}/{nfound} = {correct_dir/nfound*100:.0f}%")
pd.DataFrame(rows).to_csv("gse91061_validation_fpkm.csv", index=False)
print("Saved gse91061_validation_fpkm.csv")
