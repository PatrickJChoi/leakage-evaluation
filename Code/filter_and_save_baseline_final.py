import gzip
import pandas as pd
import numpy as np
import pickle
import time
from collections import Counter

tpm_file = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
meta_file = "GSE120575_patient_ID_single_cells.txt.gz"

start_time = time.time()

# 1. Parse metadata cells to get patient, response, therapy, etc.
meta_dict = {}
with gzip.open(meta_file, 'rt') as f:
    header = None
    for line in f:
        parts = line.strip().split('\t')
        if parts and parts[0] == "Sample name":
            header = parts
            break
    title_idx = header.index("title")
    pat_idx = [i for i, c in enumerate(header) if 'patinet' in c.lower() or 'patient' in c.lower()][0]
    resp_idx = [i for i, c in enumerate(header) if 'response' in c.lower()][0]
    therapy_idx = [i for i, c in enumerate(header) if 'therapy' in c.lower()][0]
    org_idx = header.index("organism")
    source_idx = header.index("source name")
    
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) <= max(title_idx, pat_idx, resp_idx, therapy_idx, org_idx, source_idx):
            continue
        meta_dict[parts[title_idx]] = {
            "patient": parts[pat_idx],
            "response": parts[resp_idx],
            "therapy": parts[therapy_idx],
            "organism": parts[org_idx],
            "source_name": parts[source_idx]
        }

# 2. Read TPM header to identify baseline cells
with gzip.open(tpm_file, 'rt') as f:
    tpm_barcodes = f.readline().strip().split('\t')
    tpm_patients = f.readline().strip().split('\t')

baseline_indices = [idx for idx, pat in enumerate(tpm_patients) if pat.startswith("Pre_")]
baseline_barcodes = [tpm_barcodes[idx] for idx in baseline_indices]
baseline_patients = [tpm_patients[idx] for idx in baseline_indices]

# usecols: 0 is the gene name column, and then baseline indices (offset by 1)
usecols = [0] + [idx + 1 for idx in baseline_indices]

print(f"Loading {len(baseline_indices)} baseline columns out of {len(tpm_barcodes)} from TPM file...")
df = pd.read_csv(tpm_file, sep='\t', skiprows=2, header=None, usecols=usecols)
print(f"Loaded TPM in {time.time() - start_time:.1f}s")

# Rename columns
df.columns = ['Gene'] + baseline_barcodes
df.set_index('Gene', inplace=True)

# Convert to float32
expr = df.values.astype(np.float32)
num_genes, num_baseline_cells = expr.shape
print(f"Expression matrix shape: {expr.shape}")

# Quality Checks
# 1. Mean TPM per cell (flag mean < 0.1)
cell_means = np.mean(expr, axis=0)
flagged_mean_tpm = [baseline_barcodes[i] for i in range(num_baseline_cells) if cell_means[i] < 0.1]

# 2. Genes detected per cell (flag genes < 200)
genes_detected = np.sum(expr > 0, axis=0)
flagged_low_genes = [baseline_barcodes[i] for i in range(num_baseline_cells) if genes_detected[i] < 200]

# Combine cell flags
flagged_cells = set(flagged_mean_tpm).union(set(flagged_low_genes))
print(f"Total cells flagged for low quality (mean TPM < 0.1 or genes < 200): {len(flagged_cells)}")

# 3. Identify cells belonging to patient Pre_P4
p4_cells = [baseline_barcodes[i] for i in range(num_baseline_cells) if baseline_patients[i] == "Pre_P4"]
print(f"Total cells associated with patient Pre_P4: {len(p4_cells)}")

# Filter cells: exclude flagged cells and Pre_P4 cells
cells_to_remove = flagged_cells.union(set(p4_cells))
clean_cell_indices = [i for i in range(num_baseline_cells) if baseline_barcodes[i] not in cells_to_remove]
clean_barcodes = [baseline_barcodes[i] for i in clean_cell_indices]
clean_patients = [baseline_patients[i] for i in clean_cell_indices]
num_clean_cells = len(clean_cell_indices)

print(f"New clean cell count: {num_clean_cells} (removed {len(cells_to_remove)} cells total: {len(flagged_cells)} low quality and/or {len(p4_cells)} Pre_P4)")

# Subset expression matrix to clean cells
expr_clean_cells = expr[:, clean_cell_indices]

# 4. Apply gene expression filter on the remaining cells
# Keep only genes with TPM > 1 in at least 3 remaining clean cells
surviving_genes_filter = np.sum(expr_clean_cells > 1, axis=1) >= 3
surviving_gene_names = df.index[surviving_genes_filter]
expr_final = expr_clean_cells[surviving_genes_filter, :]
num_final_genes = len(surviving_gene_names)

print(f"Final clean gene count (TPM > 1 in >= 3 clean cells): {num_final_genes}")
print(f"Final clean matrix dimensions (genes x cells): {num_final_genes} x {num_clean_cells}")

# 5. Create pandas DataFrames for final expression matrix and metadata
expr_df_clean = pd.DataFrame(
    data=expr_final,
    index=surviving_gene_names,
    columns=clean_barcodes,
    dtype=np.float32
)

metadata_records = []
for b in clean_barcodes:
    rec = {"barcode": b}
    if b in meta_dict:
        rec.update(meta_dict[b])
    else:
        # Fallback if not in metadata (though all barcodes should be)
        # We can extract patient from baseline_patients
        pat = tpm_patients[tpm_barcodes.index(b)]
        rec.update({
            "patient": pat,
            "response": "Unknown",
            "therapy": "Unknown",
            "organism": "Homo sapiens",
            "source_name": "Melanoma single cell"
        })
    metadata_records.append(rec)

metadata_df_clean = pd.DataFrame(metadata_records)
metadata_df_clean.set_index("barcode", inplace=True)

# 6. Report final patient counts, responders, non-responders
unique_patients_clean = sorted(list(set(metadata_df_clean["patient"])))
print(f"\nUnique patients remaining: {len(unique_patients_clean)}")
print("Remaining patient list:", unique_patients_clean)

patient_responses = {}
for b, row in metadata_df_clean.iterrows():
    patient_responses[row["patient"]] = row["response"]

response_counts = Counter(patient_responses.values())
print("Patient-level response counts:")
for r, c in response_counts.items():
    print(f"  - {r}: {c}")

# 7. Save to disk as pickle files
expr_df_clean.to_pickle("baseline_expression_clean.pkl")
metadata_df_clean.to_pickle("baseline_metadata_clean.pkl")
print("\nSaved expression matrix to baseline_expression_clean.pkl")
print("Saved metadata to baseline_metadata_clean.pkl")

# Verify loading works
with open("baseline_expression_clean.pkl", "rb") as f:
    verify_expr = pickle.load(f)
with open("baseline_metadata_clean.pkl", "rb") as f:
    verify_meta = pickle.load(f)
    
print(f"\nVerification:")
print(f"  - Loaded expression shape: {verify_expr.shape}")
print(f"  - Loaded metadata shape: {verify_meta.shape}")
print("Done!")
