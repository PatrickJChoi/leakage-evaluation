"""
Verified melanoma LOO-CV + nested permutation.
Mirrors ml_pipeline.py Steps 2-4 EXACTLY, but loads the already-built
pseudobulk_matrix.pkl (baseline_expression_clean.pkl is unpicklable under
this pandas version). No expected numbers referenced; honest recompute.
"""
import pandas as pd, numpy as np, pickle, time
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix

df = pd.read_pickle("pseudobulk_matrix.pkl")
gene_names = df.columns[:-1].values
X = df.iloc[:, :-1].values.astype(np.float64)
y = np.array([1 if r == "Responder" else 0 for r in df["response"]])
patient_labels = np.array(df.index)
n = len(y); k = 5
print(f"Loaded pseudobulk: X={X.shape}, labels={dict(Counter(df['response']))}")

# ---- LOO-CV (identical to ml_pipeline.py) ----
preds, probs, fold_features = [], [], []
for i in range(n):
    mask = np.ones(n, dtype=bool); mask[i] = False
    X_train, y_train = X[mask], y[mask]
    X_test = X[i:i+1]
    f_scores, _ = f_classif(X_train, y_train)
    f_scores = np.nan_to_num(f_scores, nan=-np.inf)
    top_k_idx = np.argsort(f_scores)[-k:][::-1]
    fold_features.append(gene_names[top_k_idx].tolist())
    clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    clf.fit(X_train[:, top_k_idx], y_train)
    prob = clf.predict_proba(X_test[:, top_k_idx])[0, 1]
    preds.append(int(prob >= 0.5)); probs.append(prob)
preds = np.array(preds); probs = np.array(probs)
observed_acc = accuracy_score(y, preds)
observed_auc = roc_auc_score(y, probs)
cm = confusion_matrix(y, preds, labels=[0, 1])
feature_freq = Counter([f for fold in fold_features for f in fold])
print(f"LOO acc={observed_acc:.4f} ({int(observed_acc*n)}/{n})  AUC={observed_auc:.4f}")
print("Confusion [rows=true NR,R; cols=pred NR,R]:", cm.tolist())
print("Top features:", feature_freq.most_common(8))

# ---- Permutation (identical to ml_pipeline.py Step 4) ----
rng = np.random.RandomState(42)
n_perm = 500
perm_accs = np.zeros(n_perm)
for p_idx in range(n_perm):
    y_shuf = rng.permutation(y)
    pp = np.zeros(n, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool); mask[i] = False
        X_tr, y_tr = X[mask], y_shuf[mask]
        X_te = X[i:i+1]
        f_sc, _ = f_classif(X_tr, y_tr)
        f_sc = np.nan_to_num(f_sc, nan=-np.inf)
        idx_k = np.argsort(f_sc)[-k:][::-1]
        clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        clf.fit(X_tr[:, idx_k], y_tr)
        pp[i] = int(clf.predict_proba(X_te[:, idx_k])[0, 1] >= 0.5)
    perm_accs[p_idx] = accuracy_score(y_shuf, pp)
    if (p_idx + 1) % 100 == 0:
        print(f"  perm {p_idx+1}/{n_perm} mean={np.mean(perm_accs[:p_idx+1]):.4f}")

empirical_p = np.mean(perm_accs >= observed_acc)
ge = int(np.sum(perm_accs >= observed_acc))
print(f"\nMELANOMA RESULT: acc={observed_acc:.4f} auc={observed_auc:.4f} "
      f"p={empirical_p:.4f} ({ge}/{n_perm})  mean_perm={np.mean(perm_accs):.4f} max_perm={np.max(perm_accs):.4f}")

out = {"predictions": preds, "probabilities": probs, "true_labels": y,
       "patient_labels": patient_labels, "fold_features": fold_features,
       "feature_freq": dict(feature_freq), "accuracy": observed_acc,
       "auc": observed_auc, "confusion_matrix": cm,
       "perm_accuracies": perm_accs, "p_value": empirical_p, "n_ge": ge}
with open("melanoma_verified_results_FROM_VERIFIED_SCRIPT.pkl", "wb") as f:
    pickle.dump(out, f)
print("Saved melanoma_verified_results_FROM_VERIFIED_SCRIPT.pkl")
