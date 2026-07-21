import pickle
import pandas as pd
import numpy as np
from collections import Counter

# Reload the pseudobulk matrix and re-run the LOO-CV feature selection
# to extract the exact top 20 genes programmatically

from sklearn.feature_selection import f_classif

expr = pd.read_pickle("baseline_expression_clean.pkl")
meta = pd.read_pickle("baseline_metadata_clean.pkl")

# Rebuild pseudobulk
patients = sorted(meta["patient"].unique())
gene_names = expr.index.values

expr_log = np.log2(expr.values.astype(np.float64) + 1)
expr_log_df = pd.DataFrame(expr_log, index=expr.index, columns=expr.columns)

pseudobulk_rows = []
patient_responses = []
for pat in patients:
    cells = meta.index[meta["patient"] == pat].tolist()
    cells_in_expr = [c for c in cells if c in expr_log_df.columns]
    pat_mean = expr_log_df[cells_in_expr].mean(axis=1).values
    pseudobulk_rows.append(pat_mean)
    patient_responses.append(meta.loc[cells_in_expr[0], "response"])

X = np.array(pseudobulk_rows, dtype=np.float64)
y = np.array([1 if r == "Responder" else 0 for r in patient_responses])
n = len(y)
k = 5

# Run LOO-CV feature selection only
all_features = []
for i in range(n):
    mask = np.ones(n, dtype=bool)
    mask[i] = False
    X_train, y_train = X[mask], y[mask]
    f_scores, _ = f_classif(X_train, y_train)
    f_scores = np.nan_to_num(f_scores, nan=-np.inf)
    top_k_idx = np.argsort(f_scores)[-k:][::-1]
    top_k_names = gene_names[top_k_idx]
    all_features.extend(top_k_names.tolist())

feature_freq = Counter(all_features)
print("Top 20 most frequently selected genes across 18 LOO-CV folds:")
print(f"{'Rank':<6}{'Gene':<25}{'Folds':<10}")
print("-" * 41)
for rank, (gene, count) in enumerate(feature_freq.most_common(20), 1):
    print(f"{rank:<6}{gene:<25}{count}/{n}")

# Save the top 20 gene list
top20 = [g for g, _ in feature_freq.most_common(20)]
with open("top20_genes.txt", "w") as f:
    for g in top20:
        f.write(g + "\n")
print(f"\nSaved top 20 gene list to top20_genes.txt")
