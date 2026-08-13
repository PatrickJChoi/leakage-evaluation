import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score, accuracy_score

# Load pseudobulk matrix (now in-repo, alongside this script)
df = pd.read_pickle(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Data", "pseudobulk_matrix.pkl"))
gene_names = df.columns[:-1].values
X = df.iloc[:, :-1].values.astype(np.float64)
y = np.array([1 if r == "Responder" else 0 for r in df["response"]])

n = len(y)
k = 5
predictions = []
true_labels = []
probabilities = []

for i in range(n):
    mask = np.ones(n, dtype=bool)
    mask[i] = False
    X_train, y_train = X[mask], y[mask]
    X_test = X[i:i+1]
    
    # Feature selection on training set only (ANOVA F-test)
    f_scores, _ = f_classif(X_train, y_train)
    f_scores = np.nan_to_num(f_scores, nan=-np.inf)
    top_k_idx = np.argsort(f_scores)[-k:][::-1]
    
    # Train RF on selected features
    clf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    clf.fit(X_train[:, top_k_idx], y_train)
    
    prob = clf.predict_proba(X_test[:, top_k_idx])[0, 1]
    pred = int(prob >= 0.5)
    
    predictions.append(pred)
    true_labels.append(y[i])
    probabilities.append(prob)

acc = accuracy_score(true_labels, predictions)
auc = roc_auc_score(true_labels, probabilities)
print(f"Accuracy: {acc:.4f} ({sum(np.array(predictions) == np.array(true_labels))}/{n})")
print(f"AUC: {auc:.4f}")

# Save data for plotting
np.savez("melanoma_roc_data.npz", y_true=np.array(true_labels), y_probs=np.array(probabilities))
print("Saved melanoma_roc_data.npz successfully.")
