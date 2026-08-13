"""
Computes the exact permutation p-value (b+1)/(m+1) for the melanoma LOO-CV
model directly from melanoma_verified_results.pkl, and writes the result
plus both input values to Results/melanoma_pvalue_exact.json.
"""
import json
import os
import pickle

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # papers/leakage-evaluation
MELANOMA_RESULTS_PKL = os.path.join(REPO_ROOT, "melanoma_verified_results.pkl")
OUT_PATH = os.path.join(REPO_ROOT, "Results", "melanoma_pvalue_exact.json")

with open(MELANOMA_RESULTS_PKL, "rb") as f:
    results = pickle.load(f)

n_ge = int(results["n_ge"])
perm_accuracies = results["perm_accuracies"]
m = len(perm_accuracies)

p_value_exact = (n_ge + 1) / (m + 1)

output = {
    "n_ge": n_ge,
    "m": m,
    "p_value_exact": p_value_exact,
}

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(output, f, indent=2)

print(f"n_ge = {n_ge}")
print(f"m = {m}")
print(f"p_value_exact = (n_ge + 1) / (m + 1) = {p_value_exact}")
print(f"Wrote {OUT_PATH}")
