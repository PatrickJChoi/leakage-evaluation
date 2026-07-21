"""
IBD (GSE282122) verified re-run — NumPy/parallel reimplementation of pipeline.py.
Same nested LOO-CV + genuine nested permutation logic, but memory-light so the
500-iteration permutation actually completes. Stage 1 caches the 36xG feature
matrices to npz; Stage 2 runs observed + permutation from that cache.
"""
import os, time, pickle
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings("ignore")

CACHE = "ibd_feature_cache.npz"
t0 = time.time()

if not os.path.exists(CACHE):
    import scanpy as sc
    print("Loading h5ad (one time)...")
    adata = sc.read_h5ad('../Data/GSE282122/taurus_lightweight.h5ad')
    assert adata.X.min() >= 0
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    props = adata.obs.groupby(['Patient','major'], observed=True).size().unstack(fill_value=0)
    props = props.div(props.sum(axis=1), axis=0)
    prop_cols = ['prop_'+str(c) for c in props.columns]
    patients = list(props.index)
    Xg = np.vstack([np.ravel(adata[adata.obs['Patient']==p].X.mean(axis=0)) for p in patients])
    y = adata.obs.groupby('Patient', observed=True)['Remission_status'].first().map(
        {'Remission':1,'Non_Remission':0}).reindex(patients).values.astype(int)
    genes = list(adata.var_names)
    feat_names = np.array(prop_cols + genes, dtype=object)
    Xp = props.values
    np.savez(CACHE, Xp=Xp, Xg=Xg, y=y, feat_names=feat_names, patients=np.array(patients, dtype=object))
    del adata
    print(f"Cached features ({time.time()-t0:.0f}s)")

d = np.load(CACHE, allow_pickle=True)
Xp, Xg, y, feat_names, patients = d['Xp'], d['Xg'], d['y'], d['feat_names'], d['patients']
n = len(y); n_prop = Xp.shape[1]
print(f"patients={n} prop_feats={n_prop} genes={Xg.shape[1]} labels={{0:{int((y==0).sum())},1:{int((y==1).sum())}}}")

def f_test_vec(X, yy):
    # one-way ANOVA F per column, 2 classes
    n0 = np.sum(yy==0); n1 = np.sum(yy==1); N = len(yy)
    m0 = X[yy==0].mean(0); m1 = X[yy==1].mean(0); mg = X.mean(0)
    ssb = n0*(m0-mg)**2 + n1*(m1-mg)**2
    ssw = ((X[yy==0]-m0)**2).sum(0) + ((X[yy==1]-m1)**2).sum(0)
    ssw = np.where(ssw==0, 1e-9, ssw)
    return np.nan_to_num((ssb/1)/(ssw/(N-2)), nan=-1.0)

def loo_predict(yv, k=5, track=False):
    """One full nested LOO-CV pass; returns preds (and per-fold feats if track)."""
    preds = np.zeros(n, int); probs = np.zeros(n); feats_all = []
    for i in range(n):
        tr = np.arange(n) != i
        Xg_tr = Xg[tr]; Xg_te = Xg[i:i+1]
        mu = Xg_tr.mean(0); sd = Xg_tr.std(0, ddof=1); sd[sd==0]=1.0  # ddof=1 matches pandas/pipeline.py
        Xg_tr_s = np.clip((Xg_tr-mu)/sd, None, 10)
        Xg_te_s = np.clip((Xg_te-mu)/sd, None, 10)
        Xtr = np.hstack([Xp[tr], Xg_tr_s]); Xte = np.hstack([Xp[i:i+1], Xg_te_s])
        f, _ = f_classif(Xtr, yv[tr])          # sklearn f_classif, identical to pipeline.py
        f = np.nan_to_num(f, nan=-1.0)
        idx = np.sort(np.argsort(f)[::-1][:k])  # restore ORIGINAL column order (matches SelectKBest.get_support)
        rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        rf.fit(Xtr[:, idx], yv[tr])
        preds[i] = rf.predict(Xte[:, idx])[0]
        probs[i] = rf.predict_proba(Xte[:, idx])[0,1]
        if track: feats_all.append(feat_names[idx].tolist())
    return preds, probs, feats_all

# --- observed ---
preds, probs, fold_feats = loo_predict(y, track=True)
acc = accuracy_score(y, preds); auc = roc_auc_score(y, probs)
cm = confusion_matrix(y, preds, labels=[0,1])
from collections import Counter
consensus = Counter([f for fold in fold_feats for f in fold])
print(f"LOO acc={acc:.4f} auc={auc:.4f} cm={cm.tolist()}")
print("consensus:", consensus.most_common(12))

# --- permutation (genuine nested; shuffle labels, re-run everything) ---
def one_perm(seed):
    rng = np.random.RandomState(seed)
    ys = rng.permutation(y)
    p,_,_ = loo_predict(ys)
    return float(np.mean(p==ys))

print("Running 500 nested permutations (parallel)...")
tp=time.time()
null = Parallel(n_jobs=-1, verbose=5)(delayed(one_perm)(s) for s in range(500))
null = np.array(null)
p_value = (np.sum(null>=acc)+1)/(len(null)+1)
print(f"\nIBD RESULT: acc={acc:.4f} auc={auc:.4f} p={p_value:.4f} (>=: {int(np.sum(null>=acc))}/500) "
      f"meannull={null.mean():.4f} maxnull={null.max():.4f}  ({time.time()-tp:.0f}s)")

out = {"predictions":preds,"probabilities":probs,"true_labels":y,"patients":list(patients),
       "fold_features":fold_feats,"consensus_counts":dict(consensus),"accuracy":acc,"auc":auc,
       "confusion_matrix":cm,"perm_accuracies":null,"p_value":p_value,"n_patients":n}
with open("ibd_correct_results.pkl","wb") as f:
    pickle.dump(out, f, protocol=4)
print(f"Saved ibd_correct_results.pkl (total {time.time()-t0:.0f}s)")
