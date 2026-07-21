"""
Melanoma Case Study 2 recomputed FROM RAW GSE120575 TPM text (no pickles).
Replicates filter_and_save_baseline_final.py QC + ml_pipeline.py LOO/perm.
Honest recompute; no expected numbers referenced.
"""
import gzip, time, pickle
import numpy as np
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix

tpm_file = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
meta_file = "GSE120575_patient_ID_single_cells.txt.gz"
t0 = time.time()

# --- metadata ---
meta_dict = {}
with gzip.open(meta_file, 'rt') as f:
    header = None
    for line in f:
        parts = line.strip().split('\t')
        if parts and parts[0] == "Sample name":
            header = parts; break
    title_idx = header.index("title")
    pat_idx = [i for i,c in enumerate(header) if 'patinet' in c.lower() or 'patient' in c.lower()][0]
    resp_idx = [i for i,c in enumerate(header) if 'response' in c.lower()][0]
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) <= max(title_idx, pat_idx, resp_idx): continue
        meta_dict[parts[title_idx]] = {"patient": parts[pat_idx], "response": parts[resp_idx]}

# --- TPM header: baseline (Pre_) cells ---
with gzip.open(tpm_file, 'rt') as f:
    tpm_barcodes = f.readline().strip().split('\t')
    tpm_patients = f.readline().strip().split('\t')
baseline_indices = [i for i,p in enumerate(tpm_patients) if p.startswith("Pre_")]
baseline_barcodes = [tpm_barcodes[i] for i in baseline_indices]
baseline_patients = [tpm_patients[i] for i in baseline_indices]
print(f"Baseline (Pre_) cells: {len(baseline_indices)}")

# --- load TPM for baseline cells ---
import pandas as pd
usecols = [0] + [i+1 for i in baseline_indices]
df = pd.read_csv(tpm_file, sep='\t', skiprows=2, header=None, usecols=usecols)
df.columns = ['Gene'] + baseline_barcodes
df.set_index('Gene', inplace=True)
expr = df.values.astype(np.float32)
ng, nc = expr.shape
print(f"Baseline matrix genes x cells: {expr.shape}  ({time.time()-t0:.1f}s)")

# --- QC identical to filter_and_save_baseline_final.py ---
cell_means = np.mean(expr, axis=0)
genes_detected = np.sum(expr > 0, axis=0)
flagged = set(baseline_barcodes[i] for i in range(nc) if cell_means[i] < 0.1) \
        | set(baseline_barcodes[i] for i in range(nc) if genes_detected[i] < 200)
p4_cells = set(baseline_barcodes[i] for i in range(nc) if baseline_patients[i] == "Pre_P4")
print(f"Low-quality flagged: {len(flagged)}  |  Pre_P4 cells: {len(p4_cells)}")
remove = flagged | p4_cells
keep = [i for i in range(nc) if baseline_barcodes[i] not in remove]
clean_barcodes = [baseline_barcodes[i] for i in keep]
clean_patients = [baseline_patients[i] for i in keep]
expr_clean = expr[:, keep]
gene_filter = np.sum(expr_clean > 1, axis=1) >= 3
gene_names = df.index[gene_filter].values
expr_final = expr_clean[gene_filter, :]   # genes x clean_cells
print(f"Clean cells: {len(keep)}  surviving genes: {expr_final.shape[0]}")

# --- pseudobulk: mean of log2(TPM+1) per patient ---
expr_log = np.log2(expr_final + 1)   # genes x cells
patients = sorted(set(clean_patients))
cp = np.array(clean_patients)
resp_by_pat = {}
for b, p in zip(clean_barcodes, clean_patients):
    resp_by_pat[p] = meta_dict.get(b, {}).get("response", "NA")
X = np.vstack([expr_log[:, cp == p].mean(axis=1) for p in patients])  # patients x genes
responses = [resp_by_pat[p] for p in patients]
print("Patients:", len(patients), " P4 present?", "Pre_P4" in patients)
print("Raw response values:", Counter(responses))
y = np.array([1 if r == "Responder" else 0 for r in responses])
print("Label counts (1=R,0=NR):", Counter(y.tolist()))

# --- LOO-CV ---
n = len(y); k = 5
preds, probs, fold_features = [], [], []
for i in range(n):
    m = np.ones(n, bool); m[i] = False
    Xtr, ytr, Xte = X[m], y[m], X[i:i+1]
    fs, _ = f_classif(Xtr, ytr); fs = np.nan_to_num(fs, nan=-np.inf)
    idx = np.argsort(fs)[-k:][::-1]
    fold_features.append(gene_names[idx].tolist())
    clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    clf.fit(Xtr[:, idx], ytr)
    pr = clf.predict_proba(Xte[:, idx])[0, 1]
    preds.append(int(pr >= 0.5)); probs.append(pr)
preds, probs = np.array(preds), np.array(probs)
acc = accuracy_score(y, preds); auc = roc_auc_score(y, probs)
cm = confusion_matrix(y, preds, labels=[0, 1])
freq = Counter([f for fold in fold_features for f in fold])
print(f"\nLOO acc={acc:.4f} ({int(acc*n)}/{n}) auc={auc:.4f} cm={cm.tolist()}")
print("Consensus:", freq.most_common(8))

with open("melanoma_consensus_counts.pkl", "wb") as f:
    pickle.dump(dict(freq), f)
print("Saved melanoma_consensus_counts.pkl")

# --- permutation (identical logic to ml_pipeline.py) ---
rng = np.random.RandomState(42); n_perm = 500
perm = np.zeros(n_perm)
for pi in range(n_perm):
    ys = rng.permutation(y); pp = np.zeros(n, int)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        Xtr, ytr, Xte = X[m], ys[m], X[i:i+1]
        fs, _ = f_classif(Xtr, ytr); fs = np.nan_to_num(fs, nan=-np.inf)
        idx = np.argsort(fs)[-k:][::-1]
        clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        clf.fit(Xtr[:, idx], ytr)
        pp[i] = int(clf.predict_proba(Xte[:, idx])[0, 1] >= 0.5)
    perm[pi] = accuracy_score(ys, pp)
    if (pi+1) % 100 == 0: print(f"  perm {pi+1}/{n_perm} mean={perm[:pi+1].mean():.4f}")
ge = int((perm >= acc).sum()); p_emp = float((perm >= acc).mean())
print(f"\nMELANOMA(raw): acc={acc:.4f} auc={auc:.4f} p={p_emp:.4f} ({ge}/{n_perm}) "
      f"mean_perm={perm.mean():.4f} max_perm={perm.max():.4f}")

out = {"patients": patients, "responses": responses, "y": y, "gene_count": int(expr_final.shape[0]),
       "clean_cells": len(keep), "predictions": preds, "probabilities": probs,
       "fold_features": fold_features, "feature_freq": dict(freq), "accuracy": acc,
       "auc": auc, "confusion_matrix": cm, "perm_accuracies": perm, "p_value": p_emp, "n_ge": ge}
with open("melanoma_verified_results.pkl", "wb") as f:
    pickle.dump(out, f, protocol=4)
print("Saved melanoma_verified_results.pkl")
