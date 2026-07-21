import urllib.request
import json
import time
import urllib.parse
import re
import os

gse_list = [
    "GSE16879", "GSE23597", "GSE73661", "GSE83687", 
    "GSE39587", "GSE59071", "GSE92229", "GSE117903", 
    "GSE123993", "GSE107865", "GSE159034", "GSE191328", 
    "GSE186963", "GSE111761"
]

def fetch_gse_metadata(gse):
    print(f"Fetching metadata for: {gse}")
    # Search for GDS ID using GSE accession
    url_search = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}[ACCN]&retmode=json"
    try:
        req = urllib.request.Request(url_search, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if not id_list:
                print(f"No GDS ID found for {gse}")
                return None
            
            # Fetch summary
            id_str = ",".join(id_list)
            url_sum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_str}&retmode=json"
            req_sum = urllib.request.Request(url_sum, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req_sum) as response_sum:
                data_sum = json.loads(response_sum.read().decode())
                res = data_sum.get("result", {})
                for k, v in res.items():
                    if k != 'uids' and v.get("entrytype", "").upper() == "GSE":
                        return v
        return None
    except Exception as e:
        print(f"Error fetching {gse}: {e}")
        return None

results = []
for gse in gse_list:
    meta = fetch_gse_metadata(gse)
    if meta:
        results.append(meta)
    time.sleep(0.5)

print(f"\nFetched {len(results)} records.")

# Parse and refine details manually based on literature knowledge
parsed = []
for info in results:
    gse = info.get("accession") or info.get("extid")
    title = info.get("title", "")
    summary = info.get("summary", "")
    n_samples = info.get("n_samples", 0)
    pdate = info.get("pdat", "")
    year = pdate.split("/")[0] if pdate else "Unknown"
    gpl_str = info.get("gpl", "")
    
    # Infer platform type
    # GPL570 is Affymetrix Human Genome U133 Plus 2.0 (Microarray)
    # GPL10558 is Illumina HumanHT-12 V4.0 (Microarray)
    # GPL16791 is Illumina HiSeq 2500 (RNA-Seq)
    # GPL20301 is Illumina HiSeq 4000 (RNA-Seq)
    # GPL23159 is Illumina NextSeq 500 (RNA-Seq)
    # GPL13158 is Affymetrix HT HG-U133+ PM (Microarray)
    platform = "Microarray"
    if any(g in gpl_str for g in ["16791", "20301", "23159", "11154", "15456", "18573"]):
        platform = "RNA-Seq"
    platform_desc = f"{platform} (GPL{gpl_str})"
    
    # Infer disease type
    disease = "IBD"
    if "crohn" in (title + " " + summary).lower():
        disease = "CD"
    if "colitis" in (title + " " + summary).lower():
        if disease == "CD":
            disease = "IBD" # Both Crohn's and UC
        else:
            disease = "UC"
            
    # Infer tissue type
    tissue = "Colon Biopsy"
    if "pbmc" in (title + " " + summary).lower():
        tissue = "PBMC"
    elif "blood" in (title + " " + summary).lower():
        tissue = "Whole Blood"
    elif "rectal" in (title + " " + summary).lower():
        tissue = "Rectal Biopsy"
    elif "ileal" in (title + " " + summary).lower() or "ileum" in (title + " " + summary).lower():
        tissue = "Ileal Biopsy"
        
    # Infer response labels format
    # Microarray clinical response, mucosal healing, endoscopic remission
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

# Add query mapping manually based on the previous esearch mapping:
# Query 1: 'anti-TNF response Crohn's disease microarray' -> GSE16879, GSE23597, GSE92229, GSE111761
# Query 2: 'infliximab response inflammatory bowel disease baseline' -> GSE107865, GSE73661, GSE191328, GSE186963
# Query 3: 'adalimumab ulcerative colitis RNA-seq response' -> GSE117903, GSE159034
# Query 4: 'anti-TNF non-response IBD gene expression' -> GSE111761, GSE186963, GSE83687

# Write to a file in UTF-8
output_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "geo_datasets_table.md")
with open(output_file, "w", encoding="utf-8") as out:
    out.write("# Master Table of GEO Datasets for IBD anti-TNF response\n\n")
    out.write("| GEO Accession | Title | N Samples | Disease | Tissue | Platform | Response Label Format | Year Deposited |\n")
    out.write("|---|---|---|---|---|---|---|---|\n")
    for item in parsed:
        out.write(f"| {item['GSE']} | {item['Title']} | {item['n_samples']} | {item['Disease']} | {item['Tissue']} | {item['Platform']} | {item['ResponseFormat']} | {item['Year']} |\n")

print(f"\nMarkdown table successfully saved to: {output_file}")
