import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.utils import shuffle
from sklearn.metrics import roc_auc_score
import sys
import gc

# 1. GLOBAL LEAK-FREE STEP: pseudo-bulk raw mean counts across ALL available genes
print("Loading taurus_lightweight.h5ad raw counts...")
adata = sc.read_h5ad('../Data/GSE282122/taurus_lightweight.h5ad')

# Verify that adata.X contains non-negative values
X_min_val = adata.X.min()
print(f"adata.X minimum value: {X_min_val}")
assert X_min_val >= 0, "adata.X contains negative values, expected raw counts."

# Perform cell-wise normalization and log-transformation globally (leak-free)
print("Normalizing and log-transforming raw counts...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Compute patient-level cell proportions
print("Computing patient cell proportions...")
props = adata.obs.groupby(['Patient', 'major'], observed=True).size().unstack(fill_value=0)
props = props.div(props.sum(axis=1), axis=0)
props.columns = ['prop_' + str(col) for col in props.columns]

# Compute master pseudo-bulk matrix of raw mean counts per patient across ALL available genes
print("Computing master pseudo-bulk matrix across all genes...")
X_genes = []
for p in props.index:
    p_cells = adata[adata.obs['Patient'] == p]
    p_mean = p_cells.X.mean(axis=0)
    X_genes.append(np.ravel(p_mean))

X_genes = pd.DataFrame(np.array(X_genes), index=props.index, columns=adata.var_names)

# Assert patient index alignment
assert (props.index == X_genes.index).all(), "Patient index mismatch between proportions and expressions"

# Map target labels from adata before deleting it
y_patients = adata.obs.groupby('Patient', observed=True)['Remission_status'].first().map({'Remission': 1, 'Non_Remission': 0})
y_patients = y_patients.reindex(props.index)

print(f"Master pseudo-bulk matrix dimensions: {X_genes.shape}")
print(f"Target labels distribution:\n{y_patients.value_counts()}")

# Free cell-level memory immediately
del adata
gc.collect()

# Separate features into proportions and genes
X_prop = props
X_gene = X_genes

# 2. TRUE CV LOOP (Initialize Leave-One-Out based on Patient IDs)
loo = LeaveOneOut()
final_preds = []
final_probs = []
feature_selection_counts = {}

# We pre-split because we will scale and select features strictly inside the fold
print(f"\nRunning True Nested CV on {len(y_patients)} LOO-CV folds...")
for fold_idx, (train_idx, test_idx) in enumerate(loo.split(X_prop)):
    train_patients = [y_patients.index[i] for i in train_idx]
    test_patient = y_patients.index[test_idx[0]]
    
    # Split into train (35 patients) and test (1 patient)
    X_prop_train, X_prop_test = X_prop.iloc[train_idx], X_prop.iloc[test_idx]
    X_gene_train, X_gene_test = X_gene.iloc[train_idx], X_gene.iloc[test_idx]
    y_train_fold = y_patients.iloc[train_idx]
    y_test_fold = y_patients.iloc[test_idx]
    
    # Fit scaling parameters on training genes ONLY
    train_mean = X_gene_train.mean(axis=0)
    train_std = X_gene_train.std(axis=0)
    train_std[train_std == 0] = 1.0  # Prevent division by zero
    
    # Scale and clip training genes
    X_gene_train_scaled = (X_gene_train - train_mean) / train_std
    X_gene_train_scaled = pd.DataFrame(np.clip(X_gene_train_scaled, a_min=None, a_max=10), index=X_gene_train.index, columns=X_gene_train.columns)
    
    # Scale and clip test genes using training parameters
    X_gene_test_scaled = (X_gene_test - train_mean) / train_std
    X_gene_test_scaled = pd.DataFrame(np.clip(X_gene_test_scaled, a_min=None, a_max=10), index=X_gene_test.index, columns=X_gene_test.columns)
    
    # Concatenate proportions and scaled genes
    X_train_fold = pd.concat([X_prop_train, X_gene_train_scaled], axis=1)
    X_test_fold = pd.concat([X_prop_test, X_gene_test_scaled], axis=1)
    
    # Feature selection across the ENTIRE pool of candidate genes and cell proportions
    selector = SelectKBest(f_classif, k=5)
    selector.fit(X_train_fold, y_train_fold)
    selected_features = X_train_fold.columns[selector.get_support()]
    
    # Track feature selection frequency
    for f in selected_features:
        feature_selection_counts[f] = feature_selection_counts.get(f, 0) + 1
        
    X_train_sub = X_train_fold[selected_features]
    X_test_sub = X_test_fold[selected_features]
    
    # Train Random Forest model
    model = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    model.fit(X_train_sub, y_train_fold)
    
    final_preds.append(model.predict(X_test_sub)[0])
    final_probs.append(model.predict_proba(X_test_sub)[0, 1])

# Output metrics
accuracy = np.mean(np.array(final_preds) == y_patients.values)
print("\n" + "="*40)
print(f"STRICT NESTED ACCURACY: {accuracy:.1%}")

roc_auc = roc_auc_score(y_patients.values, final_probs)
print(f"STRICT NESTED ROC-AUC: {roc_auc:.3f}")

consensus = pd.Series(feature_selection_counts).sort_values(ascending=False).head(10)
print("\nCONSENSUS SIGNATURE (Most Robust Features):")
print(consensus)
print("="*40)

# Save intermediate variables and generate Figure 3
np.savez('roc_data.npz', y_true=y_patients.values, y_probs=final_probs)
import subprocess
subprocess.run([sys.executable, "generate_figure3.py"])
print("figure3_roc_auc.png saved to Results directory.")

# 3. PERMUTATION LOOP (completely recalculating scaling and feature selection from scratch under shuffled labels)
print("\nRunning rigorous Permutation Test (500 iterations of full strict nested LOO-CV)...")

def evaluate_permutation(i, X_prop_mat, X_gene_mat, y_true_series):
    # Shuffle target labels at patient level with fixed seed per iteration
    y_shuffled_series = shuffle(y_true_series, random_state=i)
    
    preds_perm = []
    loo_perm = LeaveOneOut()
    
    for train_idx, test_idx in loo_perm.split(X_prop_mat):
        X_prop_tr, X_prop_te = X_prop_mat.iloc[train_idx], X_prop_mat.iloc[test_idx]
        X_gene_tr, X_gene_te = X_gene_mat.iloc[train_idx], X_gene_mat.iloc[test_idx]
        y_train_shuf = y_shuffled_series.iloc[train_idx]
        y_test_shuf = y_shuffled_series.iloc[test_idx]
        
        # Recalculate scaling from scratch on training genes only
        tr_mean = X_gene_tr.mean(axis=0)
        tr_std = X_gene_tr.std(axis=0)
        tr_std[tr_std == 0] = 1.0
        
        X_gene_tr_sc = (X_gene_tr - tr_mean) / tr_std
        X_gene_tr_sc = pd.DataFrame(np.clip(X_gene_tr_sc, a_min=None, a_max=10), index=X_gene_tr.index, columns=X_gene_tr.columns)
        
        X_gene_te_sc = (X_gene_te - tr_mean) / tr_std
        X_gene_te_sc = pd.DataFrame(np.clip(X_gene_te_sc, a_min=None, a_max=10), index=X_gene_te.index, columns=X_gene_te.columns)
        
        X_tr_fold = pd.concat([X_prop_tr, X_gene_tr_sc], axis=1)
        X_te_fold = pd.concat([X_prop_te, X_gene_te_sc], axis=1)
        
        # Recalculate feature selection under shuffled labels
        sel = SelectKBest(f_classif, k=5)
        sel.fit(X_tr_fold, y_train_shuf)
        selected = X_tr_fold.columns[sel.get_support()]
        
        # Train model
        rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        rf.fit(X_tr_fold[selected], y_train_shuf)
        
        # Predict
        preds_perm.append(rf.predict(X_te_fold[selected])[0])
        
    acc = np.mean(np.array(preds_perm) == y_shuffled_series.values)
    
    if (i + 1) % 50 == 0:
        print(f"Completed {i + 1}/500 permutation iterations...")
        
    return acc

null_accs = []
for i in range(500):
    null_accs.append(evaluate_permutation(i, X_prop, X_gene, y_patients))

null_accs = np.array(null_accs)
# Strict unbiased p-value formula
final_p_value = (np.sum(null_accs >= accuracy) + 1) / (len(null_accs) + 1)

print("\n" + "="*40)
print(f"FINAL ACCURACY: {accuracy:.1%}")
print(f"FINAL P-VALUE: p = {final_p_value:.4f}")
if final_p_value < 0.05:
    print("THIS IS A STATISTICALLY SIGNIFICANT DISCOVERY.")
else:
    print("THIS IS NOT STATISTICALLY SIGNIFICANT.")
print("="*40)
