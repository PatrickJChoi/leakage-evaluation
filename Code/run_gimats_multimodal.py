"""
GIMATS validation using the MULTIMODAL feature space exactly as the Methods
describes: patient-level pseudo-bulk (5 consensus genes) + cell-type proportions.
Trained on all 36 TAURUS patients, applied to the 11 GIMATS non-responders,
default 0.5 threshold. No target numbers referenced.

Honesty note: TAURUS uses a 22-type taxonomy; GIMATS uses a different one. The
mapping below is best-effort; TAURUS cell types with no GIMATS equivalent
(incl. the CD4/CD8/unconventional T split — GIMATS only has 'alpha-beta T cell')
are set to 0 for GIMATS and reported explicitly.
"""
import numpy as np, anndata as ad, pickle
from sklearn.ensemble import RandomForestClassifier

CONSENSUS = ['TNK1','IGSF8','KCNQ1','CYP4F12','SSTR1']

# TAURUS prop feature -> list of GIMATS obs['cell_type'] labels (best-effort)
PROP_MAP = {
    'prop_B': ['B cell'],
    'prop_Plasma': ['plasma cell'],
    'prop_Innate_lymphocytes': ['innate lymphoid cell'],
    'prop_Mono_macro': ['mononuclear phagocyte'],
    'prop_Endothelium': ['endothelial cell','endothelial cell of lymphatic vessel'],
    'prop_Pericyte': ['pericyte'],
    'prop_Glial': ['glial cell'],
    'prop_Mast': ['mast cell'],
    'prop_Fibroblast': ['fibroblast'],
    'prop_Ileal_epithelium': ['enterocyte of epithelium proper of ileum','ileal goblet cell'],
    # unmappable (no GIMATS equivalent / cannot split): CD4_T, CD8_T, Unconventional_T,
    # DC, Cycling_MNP, Cycling_stroma, Non_ileal_epithelium, and 5 fibroblast subtypes
}

# ---------- TAURUS training (cached CP10k+log1p pseudobulk + proportions) ----------
d = np.load('ibd_feature_cache.npz', allow_pickle=True)
feat = d['feat_names'].astype(str)
prop_cols = [f for f in feat if f.startswith('prop_')]
nprop = len(prop_cols)
gene_names = list(feat[nprop:])
Xp = d['Xp']                       # 36 x 22 proportions
Xg = d['Xg']                       # 36 x 33075 gene means
y = d['y'].astype(int)
gidx = [gene_names.index(g) for g in CONSENSUS]
Xg5 = Xg[:, gidx]                  # 36 x 5

# multimodal feature matrix = [22 proportions | 5 consensus genes]
X_tau = np.hstack([Xp, Xg5])
feat_all = prop_cols + CONSENSUS
print(f"TAURUS multimodal features: {len(prop_cols)} proportions + {len(CONSENSUS)} genes = {X_tau.shape[1]}")

# standardize on TAURUS stats, clip 10 (identical to pipeline preprocessing of features)
mu = X_tau.mean(0); sd = X_tau.std(0, ddof=1); sd[sd==0]=1.0
X_tau_s = np.clip((X_tau - mu)/sd, None, 10)
clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
clf.fit(X_tau_s, y)
imp = sorted(zip(feat_all, clf.feature_importances_), key=lambda x:-x[1])
print("Top feature importances:", [(f, round(v,3)) for f,v in imp[:8]])

# ---------- GIMATS multimodal features ----------
a = ad.read_h5ad('../Data/GSE134809/gimats_annotated.h5ad')
fn = a.var['feature_name'].astype(str).values
gcols = [int(np.where(fn==g)[0][0]) for g in CONSENSUS]
ct = a.obs['cell_type'].astype(str).values
donors = list(a.obs['donor_id'].unique())
donor_arr = a.obs['donor_id'].values
Xexpr = a[:, gcols].X
Xexpr = Xexpr.toarray() if hasattr(Xexpr,'toarray') else np.asarray(Xexpr)

mapped = [p for p in prop_cols if p in PROP_MAP]
unmapped = [p for p in prop_cols if p not in PROP_MAP]
print(f"\nGIMATS proportion mapping: {len(mapped)}/{len(prop_cols)} TAURUS props mappable; "
      f"{len(unmapped)} set to 0 (no GIMATS equivalent):")
print("  unmapped:", [p.replace('prop_','') for p in unmapped])

X_gim = np.zeros((len(donors), len(feat_all)))
ncell=[]
for i, dn in enumerate(donors):
    m = donor_arr == dn
    ncell.append(int(m.sum()))
    total = m.sum()
    ct_d = ct[m]
    # proportions (count of mapped GIMATS labels / total cells for this donor)
    for j, p in enumerate(prop_cols):
        if p in PROP_MAP:
            cnt = np.isin(ct_d, PROP_MAP[p]).sum()
            X_gim[i, j] = cnt/total
    # 5 consensus gene pseudobulk means
    X_gim[i, nprop:] = Xexpr[m].mean(0)

# fraction of GIMATS cells represented by mapped proportions (per donor)
mapped_labels = set(l for p in PROP_MAP for l in PROP_MAP[p])
frac_repr = [np.isin(ct[donor_arr==dn], list(mapped_labels)).mean() for dn in donors]
print(f"Mean fraction of GIMATS cells captured by mapped proportions: {np.mean(frac_repr):.2f}")

X_gim_s = np.clip((X_gim - mu)/sd, None, 10)
proba = clf.predict_proba(X_gim_s)[:, 1]
preds = (proba >= 0.5).astype(int)
correct = int((preds==0).sum()); sens = correct/len(donors)

print("\n--- GIMATS MULTIMODAL VALIDATION (Methods as written) ---")
print(f"{'donor':<8}{'n_cells':>8}{'P(Rem)':>9}{'pred':>9}")
for dn,nc,pb,pr in zip(donors,ncell,proba,preds):
    print(f"{dn:<8}{nc:>8}{pb:>9.3f}{('Rem' if pr else 'NonRem'):>9}")
print(f"\nPatients: {len(donors)}  correctly non-responder: {correct}/{len(donors)}  "
      f"sensitivity={sens*100:.1f}%")
print(f"Mean predicted remission probability: {proba.mean():.4f}  range [{proba.min():.3f},{proba.max():.3f}]")

out = {"method":"multimodal (22 cell-type proportions + 5 consensus genes) RF, "
       "trained on 36 TAURUS, applied to 11 GIMATS, threshold 0.5",
       "features":feat_all, "consensus_genes":CONSENSUS,
       "props_mapped":mapped, "props_unmapped":unmapped,
       "mean_frac_cells_represented":float(np.mean(frac_repr)),
       "feature_importances":dict(zip(feat_all, clf.feature_importances_.tolist())),
       "donors":donors, "n_cells":ncell, "remission_proba":proba, "preds":preds,
       "sensitivity":sens, "mean_remission_proba":float(proba.mean())}
with open("gimats_multimodal_results.pkl","wb") as f:
    pickle.dump(out, f, protocol=4)
print("Saved gimats_multimodal_results.pkl")
