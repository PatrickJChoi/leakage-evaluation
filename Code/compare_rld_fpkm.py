import gzip, pandas as pd, numpy as np
from collections import Counter
from scipy.stats import mannwhitneyu

reported = pd.read_csv("gse91061_validation_results.csv")

def load_meta():
    metadata = {}
    with gzip.open("GSE91061_series_matrix.txt.gz",'rt',encoding='utf-8',errors='ignore') as f:
        for line in f:
            ls=line.strip()
            if ls.startswith('!'):
                p=ls.split('\t'); metadata.setdefault(p[0],[]).append([v.strip('"') for v in p[1:]])
            elif ls.startswith('"ID_REF"'): break
    titles=metadata["!Sample_title"][0]; chars=metadata.get("!Sample_characteristics_ch1",[])
    recs=[]
    for i in range(len(titles)):
        r={"title":titles[i]}
        for cl in chars:
            if i<len(cl) and ':' in cl[i]:
                r[cl[i].split(':')[0].strip().lower()]=cl[i].split(':',1)[1].strip()
        recs.append(r)
    m=pd.DataFrame(recs)
    pre=m[m['visit (pre or on treatment)']=='Pre']
    b=pre[pre['response'].isin(['PD','PRCR'])].copy()
    b['grp']=b['response'].map({'PRCR':'Responder','PD':'Non-responder'})
    return b

def pvals_from(path):
    expr=pd.read_csv(path,index_col=0,compression='gzip'); expr.index=[int(x) for x in expr.index]
    b=load_meta(); avail=[s for s in b['title'] if s in expr.columns]
    b=b[b['title'].isin(avail)].set_index('title')
    resp=b.index[b['grp']=='Responder'].tolist(); nonr=b.index[b['grp']=='Non-responder'].tolist()
    out={}
    for _,row in reported.iterrows():
        eid=int(row['entrez_id'])
        if eid not in expr.index: out[row['gene']]=(None,None,None); continue
        rv=expr.loc[eid,resp].values.astype(float); nv=expr.loc[eid,nonr].values.astype(float)
        U,p=mannwhitneyu(rv,nv,alternative='two-sided')
        out[row['gene']]=(p,float(np.mean(rv)),float(U))
    return out

rld=pvals_from("GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz")
fpkm=pvals_from("GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")

print(f"{'gene':<10}{'reported_p':>12}{'rld_p':>12}{'fpkm_p':>12}   matches")
n_rld=n_fpkm=n_neither=0
for _,row in reported.iterrows():
    g=row['gene']; rep=row['p_value']
    rp=rld[g][0]; fp=fpkm[g][0]
    m_rld = rp is not None and abs(rp-rep)<1e-9
    m_fpkm= fp is not None and abs(fp-rep)<1e-9
    tag = 'rld' if m_rld and not m_fpkm else ('fpkm' if m_fpkm and not m_rld else ('BOTH' if m_rld and m_fpkm else 'NEITHER'))
    if m_rld and not m_fpkm: n_rld+=1
    elif m_fpkm and not m_rld: n_fpkm+=1
    elif not m_rld and not m_fpkm: n_neither+=1
    rps = f"{rp:.6f}" if rp is not None else "NA"
    fps = f"{fp:.6f}" if fp is not None else "NA"
    print(f"{g:<10}{rep:>12.6f}{rps:>12}{fps:>12}   {tag}")
print(f"\nreported matches rld only: {n_rld}/13 | fpkm only: {n_fpkm}/13 | neither: {n_neither}/13")
