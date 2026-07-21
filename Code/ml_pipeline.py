import pandas as pd
import numpy as np
import pickle
import time
from collections import Counter

from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix

# ============================================================
# Step 1 — Pseudo-bulk aggregation
# ============================================================
print("=" * 60)
print("STEP 1: Pseudo-bulk aggregation")
print("=" * 60)

t0 = time.time()

expr = pd.read_pickle("baseline_expression_clean.pkl")   # genes × cells
meta = pd.read_pickle("baseline_metadata_clean.pkl")     # cells × metadata

print(f"Loaded expression: {expr.shape}  metadata: {meta.shape}")

# log2(TPM + 1) transform
expr_log = np.log2(expr.values.astype(np.float64) + 1)
expr_log_df = pd.DataFrame(expr_log, index=expr.index, columns=expr.columns)

# Aggregate per patient: mean across cells
patients = sorted(meta["patient"].unique())
gene_names = expr.index.values

pseudobulk_rows = []
patient_labels = []
patient_responses = []

for pat in patients:
    cells = meta.index[meta["patient"] == pat].tolist()
    cells_in_expr = [c for c in cells if c in expr_log_df.columns]
    pat_mean = expr_log_df[cells_in_expr].mean(axis=1).values
    pseudobulk_rows.append(pat_mean)
    patient_labels.append(pat)
    patient_responses.append(meta.loc[cells_in_expr[0], "response"])

X = np.array(pseudobulk_rows, dtype=np.float64)          # 18 × genes
y = np.array([1 if r == "Responder" else 0 for r in patient_responses])
patient_labels = np.array(patient_labels)

pseudobulk_df = pd.DataFrame(X, index=patient_labels, columns=gene_names)
pseudobulk_df["response"] = patient_responses
pseudobulk_df.to_pickle("pseudobulk_matrix.pkl")

print(f"Pseudo-bulk matrix shape (patients × genes): {X.shape}")
print(f"Labels: {dict(Counter(patient_responses))}")
print(f"Saved pseudobulk_matrix.pkl  ({time.time()-t0:.1f}s)")

# ============================================================
# Step 2 — Nested LOO-CV pipeline
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Nested Leave-One-Out Cross-Validation")
print("=" * 60)

n = len(y)
k = 5
predictions = []
true_labels = []
probabilities = []
fold_features = []

for i in range(n):
    # Hold out patient i
    mask = np.ones(n, dtype=bool)
    mask[i] = False
    X_train, y_train = X[mask], y[mask]
    X_test = X[i:i+1]

    # Feature selection on training set only (ANOVA F-test)
    f_scores, p_vals = f_classif(X_train, y_train)
    # Handle NaN f-scores (constant features) by setting to -inf
    f_scores = np.nan_to_num(f_scores, nan=-np.inf)
    top_k_idx = np.argsort(f_scores)[-k:][::-1]
    top_k_names = gene_names[top_k_idx]
    fold_features.append(top_k_names.tolist())

    # Train RF on selected features
    clf = RandomForestClassifier(
        n_estimators=100, max_depth=3, random_state=42
    )
    clf.fit(X_train[:, top_k_idx], y_train)

    # Predict
    prob = clf.predict_proba(X_test[:, top_k_idx])[0, 1]
    pred = int(prob >= 0.5)

    predictions.append(pred)
    true_labels.append(y[i])
    probabilities.append(prob)

    print(f"  Fold {i+1:2d}/{n}  patient={patient_labels[i]:<10}  "
          f"true={y[i]}  pred={pred}  prob={prob:.3f}  "
          f"features={list(top_k_names)}")

predictions = np.array(predictions)
true_labels = np.array(true_labels)
probabilities = np.array(probabilities)

# ============================================================
# Step 3 — Compute metrics
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Metrics")
print("=" * 60)

observed_acc = accuracy_score(true_labels, predictions)
observed_auc = roc_auc_score(true_labels, probabilities)

print(f"LOO-CV Accuracy : {observed_acc:.4f}  ({int(observed_acc*n)}/{n})")
print(f"ROC-AUC         : {observed_auc:.4f}")

# Feature selection frequency
all_features = [f for fold in fold_features for f in fold]
feature_freq = Counter(all_features)
print(f"\nFeature selection frequency (top features across {n} folds):")
for gene, count in feature_freq.most_common(20):
    print(f"  {gene:<20s}  selected in {count}/{n} folds")

# Confusion matrix
cm = confusion_matrix(true_labels, predictions, labels=[0, 1])
print(f"\nConfusion Matrix (rows=true, cols=pred):")
print(f"                  Pred NR   Pred R")
print(f"  True NR          {cm[0,0]:5d}    {cm[0,1]:5d}")
print(f"  True R           {cm[1,0]:5d}    {cm[1,1]:5d}")

# ============================================================
# Step 4 — Permutation test  (500 permutations)
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Permutation test (500 permutations)")
print("=" * 60)

rng = np.random.RandomState(42)
n_perm = 500
perm_accs = np.zeros(n_perm)

t_perm = time.time()
for p_idx in range(n_perm):
    y_shuf = rng.permutation(y)
    perm_preds = np.zeros(n, dtype=int)

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        X_tr, y_tr = X[mask], y_shuf[mask]
        X_te = X[i:i+1]

        f_sc, _ = f_classif(X_tr, y_tr)
        f_sc = np.nan_to_num(f_sc, nan=-np.inf)
        idx_k = np.argsort(f_sc)[-k:][::-1]

        clf = RandomForestClassifier(
            n_estimators=100, max_depth=3, random_state=42
        )
        clf.fit(X_tr[:, idx_k], y_tr)
        prob_p = clf.predict_proba(X_te[:, idx_k])[0, 1]
        perm_preds[i] = int(prob_p >= 0.5)

    perm_accs[p_idx] = accuracy_score(y_shuf, perm_preds)

    if (p_idx + 1) % 50 == 0:
        elapsed = time.time() - t_perm
        print(f"  Permutation {p_idx+1:4d}/{n_perm}  "
              f"elapsed={elapsed:.1f}s  "
              f"mean_perm_acc={np.mean(perm_accs[:p_idx+1]):.4f}")

empirical_p = np.mean(perm_accs >= observed_acc)
print(f"\nPermutation test complete.")
print(f"  Observed accuracy : {observed_acc:.4f}")
print(f"  Mean perm accuracy: {np.mean(perm_accs):.4f}")
print(f"  Max  perm accuracy: {np.max(perm_accs):.4f}")
print(f"  Empirical p-value : {empirical_p:.4f}  "
      f"({int(np.sum(perm_accs >= observed_acc))}/{n_perm})")

# ============================================================
# Step 5 — Final Report
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: FINAL REPORT")
print("=" * 60)
print(f"  LOO-CV Accuracy     : {observed_acc:.4f}  ({int(observed_acc*n)}/{n})")
print(f"  ROC-AUC             : {observed_auc:.4f}")
print(f"  Permutation p-value : {empirical_p:.4f}  (500 permutations)")
print()
print(f"  Confusion Matrix:")
print(f"                       Pred NR   Pred R")
print(f"    True NR             {cm[0,0]:5d}    {cm[0,1]:5d}")
print(f"    True R              {cm[1,0]:5d}    {cm[1,1]:5d}")
print()
print(f"  Top features by selection frequency:")
for rank, (gene, count) in enumerate(feature_freq.most_common(20), 1):
    print(f"    {rank:2d}. {gene:<20s}  {count}/{n} folds")
print()
print(f"  Per-patient predictions:")
for i in range(n):
    status = "OK" if predictions[i] == true_labels[i] else "MISS"
    label_str = "R" if true_labels[i] == 1 else "NR"
    pred_str = "R" if predictions[i] == 1 else "NR"
    print(f"    {patient_labels[i]:<10}  true={label_str:>2}  "
          f"pred={pred_str:>2}  prob={probabilities[i]:.3f}  {status}")

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
