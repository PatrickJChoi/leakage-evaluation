"""
Honest GIMATS (GSE134809) external validation, consistent with the manuscript's
STATED method (Leakage_Evaluation_Preprint, Case Study 1 External Validation):
  - consensus 5-feature (5-gene) model
  - RandomForest(100 trees, max_depth=3, random_state=42)
  - trained on the FULL TAURUS discovery cohort (all 36 patients)
  - applied to the 11 GIMATS non-responders (pseudobulk, identical preprocessing)
Consensus genes = top-5 by LOO selection frequency recomputed from the pipeline:
  TNK1, IGSF8, KCNQ1, CYP4F12, SSTR1.
No target numbers referenced; reports whatever it gets.
"""
import numpy as np, anndata as ad, pickle
from sklearn.ensemble import RandomForestClassifier

CONSENSUS = ['TNK1','IGSF8','KCNQ1','CYP4F12','SSTR1']

# ---------- TAURUS training data (cached: CP10k+log1p pseudobulk means) ----------
d = np.load('ibd_feature_cache.npz', allow_pickle=True)
feat = d['feat_names'].astype(str)
nprop = int(np.char.startswith(feat, 'prop_').sum())
gene_names = feat[nprop:]
Xg = d['Xg']                       # 36 x 33075  (patient pseudobulk, log1p CP10k)
y = d['y'].astype(int)             # 1=Remission, 0=Non_Remission
gidx = [list(gene_names).index(g) for g in CONSENSUS]
X_tau = Xg[:, gidx]                # 36 x 5
print(f"TAURUS: {X_tau.shape[0]} patients, genes={CONSENSUS}, labels R/NR = {int((y==1).sum())}/{int((y==0).sum())}")

# standardize on TAURUS (train) stats, clip at 10 — identical to pipeline.py
mu = X_tau.mean(0); sd = X_tau.std(0, ddof=1); sd[sd==0]=1.0
X_tau_s = np.clip((X_tau - mu)/sd, None, 10)
clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
clf.fit(X_tau_s, y)   # single final model on all 36 patients

# ---------- GIMATS: pseudobulk the 11 non-responders ----------
a = ad.read_h5ad('../Data/GSE134809/gimats_annotated.h5ad')  # in memory
fn = a.var['feature_name'].astype(str).values
col_for = {}
for g in CONSENSUS:
    w = np.where(fn == g)[0]
    col_for[g] = int(w[0]) if len(w) else None
    if col_for[g] is None: print("WARNING: gene missing in GIMATS:", g)
cols = [col_for[g] for g in CONSENSUS]

donors = list(a.obs['donor_id'].unique())
X_gim = np.zeros((len(donors), len(CONSENSUS)))
Xmat = a[:, cols].X
Xmat = Xmat.toarray() if hasattr(Xmat, 'toarray') else np.asarray(Xmat)
donor_arr = a.obs['donor_id'].values
ncell = []
for i, dn in enumerate(donors):
    m = donor_arr == dn
    ncell.append(int(m.sum()))
    X_gim[i] = Xmat[m].mean(0)      # pseudobulk mean of log1p CP10k, identical to TAURUS
print(f"GIMATS: {len(donors)} donors (all non-responders by cohort definition), cells={ncell}")

# apply TAURUS standardization, predict
X_gim_s = np.clip((X_gim - mu)/sd, None, 10)
proba = clf.predict_proba(X_gim_s)[:, 1]         # P(Remission)
preds = (proba >= 0.5).astype(int)               # paper's decision threshold (default 0.5)
# ground truth: all GIMATS = Non_Remission (0). Sensitivity = correctly called non-responder.
correct_nonresp = int((preds == 0).sum())
sensitivity = correct_nonresp / len(donors)

print("\n--- GIMATS EXTERNAL VALIDATION (paper's stated method) ---")
print(f"{'donor':<8}{'n_cells':>8}{'P(Remission)':>14}{'pred':>8}")
for dn, nc, pb, pr in zip(donors, ncell, proba, preds):
    print(f"{dn:<8}{nc:>8}{pb:>14.3f}{('Rem' if pr else 'NonRem'):>8}")
print(f"\nPatients tested: {len(donors)}")
print(f"Correctly classified as non-responders (pred=NonRem): {correct_nonresp}/{len(donors)}  -> sensitivity = {sensitivity*100:.1f}%")
print(f"Mean predicted remission probability: {proba.mean():.4f}")
print(f"Range of remission probability: [{proba.min():.3f}, {proba.max():.3f}]")

out = {"method":"consensus 5-gene RF, trained on full TAURUS (36 pts), applied to GIMATS 11",
       "consensus_genes":CONSENSUS, "donors":donors, "n_cells":ncell,
       "remission_proba":proba, "preds":preds, "sensitivity":sensitivity,
       "mean_remission_proba":float(proba.mean()), "tau_mu":mu, "tau_sd":sd}
with open("gimats_correct_results.pkl","wb") as f:
    pickle.dump(out, f, protocol=4)
print("Saved gimats_correct_results.pkl")
