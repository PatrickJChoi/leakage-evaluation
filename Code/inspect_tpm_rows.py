import gzip

tpm_file = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"

with gzip.open(tpm_file, 'rt') as f:
    for i in range(10):
        line = f.readline().strip()
        parts = line.split('\t')
        print(f"Row {i+1} (first 6 columns): {parts[:6]}")
