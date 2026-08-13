import urllib.request
import json
import time
import urllib.parse
import os

gse_list = [
    "GSE16879", "GSE23597", "GSE73661", "GSE107865", 
    "GSE159034", "GSE191328", "GSE186963", "GSE111761", 
    "GSE14580", "GSE12251", "GSE59071"
]

def fetch_gse_exact(gse):
    print(f"Fetching metadata for: {gse}")
    # Search in GDS database for this GSE
    url_search = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}[ACCN]&retmode=json"
    try:
        req = urllib.request.Request(url_search, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if not id_list:
                return None
            
            # Fetch summary
            id_str = ",".join(id_list)
            url_sum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_str}&retmode=json"
            req_sum = urllib.request.Request(url_sum, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req_sum) as response_sum:
                data_sum = json.loads(response_sum.read().decode())
                res = data_sum.get("result", {})
                
                # Look for the exact GSE accession matching entry
                for k, v in res.items():
                    if k == 'uids':
                        continue
                    # Match accession exactly and make sure entrytype is GSE
                    acc = v.get("accession") or v.get("extid")
                    entrytype = v.get("entrytype", "").upper()
                    if acc == gse and entrytype == "GSE":
                        return v
        return None
    except Exception as e:
        print(f"Error fetching {gse}: {e}")
        return None

results = []
for gse in gse_list:
    meta = fetch_gse_exact(gse)
    if meta:
        results.append(meta)
    else:
        print(f"Failed to find exact GSE for {gse}")
    time.sleep(0.5)

print(f"\nFetched {len(results)} records.")

parsed = []
for info in results:
    gse = info.get("accession") or info.get("extid")
    title = info.get("title", "")
    summary = info.get("summary", "")
    n_samples = info.get("n_samples", 0)
    pdate = info.get("pdat", "")
    year = pdate.split("/")[0] if pdate else "Unknown"
    gpl_str = info.get("gpl", "")
    
    # Platform type
    platform = "Microarray"
    if any(g in gpl_str for g in ["16791", "20301", "23159", "11154", "15456", "18573"]):
        platform = "RNA-Seq"
    platform_desc = f"{platform} (GPL{gpl_str})"
    
    # Disease
    disease = "IBD"
    if "crohn" in (title + " " + summary).lower():
        disease = "CD"
    if "colitis" in (title + " " + summary).lower():
        if disease == "CD":
            disease = "IBD"
        else:
            disease = "UC"
            
    # Tissue
    tissue = "Colon Biopsy"
    if "pbmc" in (title + " " + summary).lower():
        tissue = "PBMC"
    elif "blood" in (title + " " + summary).lower():
        tissue = "Whole Blood"
    elif "rectal" in (title + " " + summary).lower():
        tissue = "Rectal Biopsy"
    elif "ileal" in (title + " " + summary).lower() or "ileum" in (title + " " + summary).lower():
        tissue = "Ileal Biopsy"
        
    # Response Labels
    response_format = "Response / Non-response"
    title_summary = (title + " " + summary).lower()
    if "mucosal healing" in title_summary or "healing" in title_summary:
        response_format = "Mucosal Healing"
    elif "remission" in title_summary:
        response_format = "Clinical/Endoscopic Remission"
    elif "mayo" in title_summary:
        response_format = "Mayo Score change"
        
    parsed.append({
        "GSE": gse,
        "Title": title,
        "n_samples": n_samples,
        "Disease": disease,
        "Tissue": tissue,
        "Platform": platform_desc,
        "ResponseFormat": response_format,
        "Year": year
    })

# Write to a file in UTF-8
output_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "geo_datasets_table_correct.md")
with open(output_file, "w", encoding="utf-8") as out:
    out.write("# Master Table of GEO Datasets for IBD anti-TNF response\n\n")
    out.write("| GEO Accession | Title | N Samples | Disease | Tissue | Platform | Response Label Format | Year Deposited |\n")
    out.write("|---|---|---|---|---|---|---|---|\n")
    for item in parsed:
        out.write(f"| {item['GSE']} | {item['Title']} | {item['n_samples']} | {item['Disease']} | {item['Tissue']} | {item['Platform']} | {item['ResponseFormat']} | {item['Year']} |\n")

print(f"\nMarkdown table successfully saved to: {output_file}")
