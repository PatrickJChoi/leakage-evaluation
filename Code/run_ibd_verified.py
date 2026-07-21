"""
Verified end-to-end re-run of the IBD single-cell analysis (GSE282122).
Identical logic to pipeline.py (nested LOO-CV + genuine 500-iter nested
permutation), but saves ALL outputs to a pickle and skips figure subprocess.
No expected numbers referenced; honest recompute.
"""
import scanpy as sc, pandas as pd, numpy as np, pickle, time, gc
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.utils import shuffle
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score

t0 = time.time()
print("Loading taurus_lightweight.h5ad ...")
adata = sc.read_h5ad('../Data/GSE282122/taurus_lightweight.h5ad')
assert adata.X.min() >= 0, "negative values in X"
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

props = adata.obs.groupby(['Patient','major'], observed=True).size().unstack(fill_value=0)
props = props.div(props.sum(axis=1), axis=0)
props.columns = ['prop_' + str(c) for c in props.columns]

X_genes = []
for p in props.index:
    pc = adata[adata.obs['Patient'] == p]
    X_genes.append(np.ravel(pc.X.mean(axis=0)))
X_genes = pd.DataFrame(np.array(X_genes), index=props.index, columns=adata.var_names)
assert (props.index == X_genes.index).all()

y_patients = adata.obs.groupby('Patient', observed=True)['Remission_status'].first().map(
    {'Remission':1,'Non_Remission':0}).reindex(props.index)
del adata; gc.collect()
X_prop, X_gene = props, X_genes
print(f"patients={len(y_patients)} genes={X_gene.shape[1]} labels={y_patients.value_counts().to_dict()} ({time.time()-t0:.0f}s)")

def one_fold(Xp_tr, Xp_te, Xg_tr, Xg_te, y_tr, k=5):
    m = Xg_tr.mean(axis=0); s = Xg_tr.std(axis=0); s[s==0]=1.0
    Xg_tr_s = pd.DataFrame(np.clip((Xg_tr-m)/s, None, 10), index=Xg_tr.index, columns=Xg_tr.columns)
    Xg_te_s = pd.DataFrame(np.clip((Xg_te-m)/s, None, 10), index=Xg_te.index, columns=Xg_te.columns)
    Xtr = pd.concat([Xp_tr, Xg_tr_s], axis=1); Xte = pd.concat([Xp_te, Xg_te_s], axis=1)
    sel = SelectKBest(f_classif, k=k).fit(Xtr, y_tr)
    feats = Xtr.columns[sel.get_support()]
    rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42).fit(Xtr[feats], y_tr)
    return rf.predict(Xte[feats])[0], rf.predict_proba(Xte[feats])[0,1], list(feats)

# --- observed LOO-CV ---
loo = LeaveOneOut(); preds=[]; probs=[]; fold_feats=[]; counts={}
for tr, te in loo.split(X_prop):
    pr, pb, feats = one_fold(X_prop.iloc[tr], X_prop.iloc[te], X_gene.iloc[tr], X_gene.iloc[te], y_patients.iloc[tr])
    preds.append(pr); probs.append(pb); fold_feats.append(feats)
    for f in feats: counts[f]=counts.get(f,0)+1
preds=np.array(preds); probs=np.array(probs)
acc = accuracy_score(y_patients.values, preds); auc = roc_auc_score(y_patients.values, probs)
cm = confusion_matrix(y_patients.values, preds, labels=[0,1])
consensus = pd.Series(counts).sort_values(ascending=False)
print(f"\nLOO acc={acc:.4f} auc={auc:.4f} cm={cm.tolist()}")
print("Top consensus:\n", consensus.head(12).to_dict())

# --- genuine nested permutation (identical to pipeline.py) ---
def perm_acc(i):
    ys = shuffle(y_patients, random_state=i); pp=[]
    for tr, te in LeaveOneOut().split(X_prop):
        pr,_,_ = one_fold(X_prop.iloc[tr], X_prop.iloc[te], X_gene.iloc[tr], X_gene.iloc[te], ys.iloc[tr])
        pp.append(pr)
    return np.mean(np.array(pp)==ys.values)
print("\nRunning 500 nested permutations ...")
null=[]
tp=time.time()
for i in range(500):
    null.append(perm_acc(i))
    if (i+1)%50==0: print(f"  {i+1}/500  meannull={np.mean(null):.4f}  ({time.time()-tp:.0f}s)")
null=np.array(null)
p_value = (np.sum(null>=acc)+1)/(len(null)+1)
print(f"\nIBD RESULT: acc={acc:.4f} auc={auc:.4f} p={p_value:.4f} (>= : {int(np.sum(null>=acc))}/500) meannull={null.mean():.4f} maxnull={null.max():.4f}")

out = {"predictions":preds,"probabilities":probs,"true_labels":y_patients.values,
       "patient_index":list(y_patients.index),"fold_features":fold_feats,
       "consensus_counts":dict(consensus),"accuracy":acc,"auc":auc,
       "confusion_matrix":cm,"perm_accuracies":null,"p_value":p_value,"n_patients":len(y_patients)}
with open("ibd_verified_results.pkl","wb") as f:
    pickle.dump(out, f, protocol=4)
print(f"Saved ibd_verified_results.pkl (total {time.time()-t0:.0f}s)")
