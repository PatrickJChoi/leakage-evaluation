import urllib.request
import gzip
import os

# Download GSE91061 series matrix
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz"
outfile = "GSE91061_series_matrix.txt.gz"

if not os.path.exists(outfile):
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, outfile)
    print(f"Downloaded: {os.path.getsize(outfile)} bytes")
else:
    print(f"File already exists: {outfile} ({os.path.getsize(outfile)} bytes)")

# Parse the series matrix
print("\n" + "=" * 60)
print("PARSING GSE91061 SERIES MATRIX")
print("=" * 60)

metadata = {}
data_started = False
with gzip.open(outfile, 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line_stripped = line.strip()
        if not line_stripped:
            continue
        if line_stripped.startswith('!'):
            parts = line_stripped.split('\t')
            key = parts[0]
            vals = [v.strip('"') for v in parts[1:]]
            if key not in metadata:
                metadata[key] = []
            metadata[key].append(vals)
        elif line_stripped.startswith('"ID_REF"'):
            data_started = True
            break

# Basic info
if "!Series_title" in metadata:
    print(f"Title: {metadata['!Series_title'][0][0]}")
if "!Series_summary" in metadata:
    for s in metadata["!Series_summary"]:
        print(f"Summary: {s[0][:200]}")
if "!Series_overall_design" in metadata:
    print(f"Design: {metadata['!Series_overall_design'][0][0][:200]}")
if "!Series_type" in metadata:
    print(f"Type: {metadata['!Series_type'][0][0]}")

# Sample count
geo_ids = metadata.get("!Sample_geo_accession", [[]])[0]
titles = metadata.get("!Sample_title", [[]])[0]
print(f"\nTotal samples: {len(geo_ids)}")

# Characteristics
chars = metadata.get("!Sample_characteristics_ch1", [])
print(f"\nNumber of characteristics lines: {len(chars)}")

# Print all unique characteristic keys
char_keys = set()
for char_line in chars:
    for val in char_line:
        if ':' in val:
            char_keys.add(val.split(':')[0].strip())

print(f"Characteristic keys found: {sorted(char_keys)}")

# Extract per-sample metadata
num_samples = len(geo_ids)
sample_data = []
for i in range(num_samples):
    rec = {"geo": geo_ids[i], "title": titles[i] if i < len(titles) else ""}
    for char_line in chars:
        if i < len(char_line):
            val = char_line[i]
            if ':' in val:
                k = val.split(':')[0].strip().lower()
                v = val.split(':', 1)[1].strip()
                rec[k] = v
    sample_data.append(rec)

# Print first 5 samples
print("\nFirst 5 samples:")
for s in sample_data[:5]:
    print(f"  {s}")

# Check response labels
response_key = None
for k in ["response", "best response", "recist", "clinical response", 
           "best overall response", "bor", "clinical benefit"]:
    if any(k in s for s in [list(sd.keys()) for sd in sample_data][0]):
        response_key = k
        break

# Search all keys for response-related ones
all_keys = set()
for sd in sample_data:
    all_keys.update(sd.keys())
print(f"\nAll metadata keys: {sorted(all_keys)}")

# Find response-related keys
resp_keys = [k for k in all_keys if any(r in k.lower() for r in 
             ["response", "recist", "benefit", "bor", "outcome"])]
print(f"Response-related keys: {resp_keys}")

# Print unique values for each response-related key
for rk in resp_keys:
    vals = [sd.get(rk, "MISSING") for sd in sample_data]
    from collections import Counter
    print(f"\n  Key '{rk}' unique values: {dict(Counter(vals))}")

# Check for timepoint/pre-treatment
time_keys = [k for k in all_keys if any(t in k.lower() for t in 
             ["time", "pre", "post", "baseline", "treatment", "on-treatment",
              "biopsy", "sample"])]
print(f"\nTimepoint-related keys: {time_keys}")
for tk in time_keys:
    vals = [sd.get(tk, "MISSING") for sd in sample_data]
    from collections import Counter
    print(f"  Key '{tk}' unique values: {dict(Counter(vals))}")

# Platform
if "!Sample_platform_id" in metadata:
    platforms = metadata["!Sample_platform_id"][0]
    from collections import Counter
    print(f"\nPlatforms: {dict(Counter(platforms))}")

# Supplementary files
if "!Series_supplementary_file" in metadata:
    print("\nSupplementary files:")
    for sf in metadata["!Series_supplementary_file"]:
        for f in sf:
            print(f"  {f}")
