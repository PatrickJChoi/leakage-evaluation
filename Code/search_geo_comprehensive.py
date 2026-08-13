import urllib.request
import json
import time
import urllib.parse
import re

queries = [
    "anti-TNF response Crohn's disease microarray",
    "infliximab response inflammatory bowel disease baseline",
    "adalimumab ulcerative colitis RNA-seq response",
    "anti-TNF non-response IBD gene expression"
]

known_gse = ["GSE107865", "GSE159034"]

def search_geo(query):
    print(f"Searching: {query}")
    encoded_query = urllib.parse.quote_plus(query)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={encoded_query}&retmode=json&retmax=100"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            id_list = data.get("esearchresult", {}).get("idlist", [])
            return id_list
    except Exception as e:
        print(f"Error searching Entrez: {e}")
        return []

def get_summary(ids):
    if not ids:
        return {}
    # Fetch in batches of 100 to avoid long URLs
    batch_size = 100
    results = {}
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        id_str = ",".join(batch)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_str}&retmode=json"
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode())
                batch_res = data.get("result", {})
                for k, v in batch_res.items():
                    if k != 'uids':
                        results[k] = v
            time.sleep(0.5)
        except Exception as e:
            print(f"Error fetching summaries in batch: {e}")
    return results

# Map GDS IDs to their summaries
all_results = {}
query_mappings = {}  # GSE -> list of query terms that found it

# Run search queries
for idx, q in enumerate(queries, 1):
    ids = search_geo(q)
    if ids:
        summaries = get_summary(ids)
        for uid, info in summaries.items():
            if info.get("entrytype", "").upper() == "GSE":
                gse = info.get("accession") or info.get("extid")
                if gse:
                    all_results[gse] = info
                    if gse not in query_mappings:
                        query_mappings[gse] = []
                    query_mappings[gse].append(f"Query {idx}")
    time.sleep(1)

# Add known GSEs if not already found (we will search them directly)
for gse in known_gse:
    if gse not in all_results:
        print(f"Searching known GSE directly: {gse}")
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}[ACCN]&retmode=json"
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode())
                id_list = data.get("esearchresult", {}).get("idlist", [])
                if id_list:
                    summaries = get_summary(id_list)
                    for uid, info in summaries.items():
                        all_results[gse] = info
                        query_mappings[gse] = ["Known GSE"]
        except Exception as e:
            print(f"Error searching known GSE: {e}")

# Parse summaries
parsed_data = []
for gse, info in all_results.items():
    title = info.get("title", "")
    summary = info.get("summary", "")
    n_samples = info.get("n_samples", 0)
    pdate = info.get("pdat", "")
    year = pdate.split("/")[0] if pdate else "Unknown"
    gpl = info.get("gpl", "")
    
    # We want to check if it matches the criteria:
    # 1. Studies anti-TNF response in IBD (Crohn's, UC, or IBD).
    # Keywords in title or summary: anti-TNF, infliximab, adalimumab, Remicade, Humira, Crohn, colitis, IBD, response, responder.
    title_summary = (title + " " + summary).lower()
    
    is_ibd = any(kw in title_summary for kw in ["crohn", "colitis", "ibd", "bowel disease"])
    is_antitnf = any(kw in title_summary for kw in ["anti-tnf", "infliximab", "adalimumab", "remicade", "humira", "anti-tumor necrosis factor", "golimumab", "certolizumab"])
    is_response = any(kw in title_summary for kw in ["respon", "nonrespon", "remission", "healing", "efficacy"])
    
    # We will keep all that are related, and then manual/semi-manual filtering
    # But let's check:
    is_relevant = is_ibd and is_antitnf and is_response
    
    # Include known ones regardless
    if gse in known_gse:
        is_relevant = True
        
    if is_relevant:
        parsed_data.append({
            "GSE": gse,
            "Title": title,
            "n_samples": n_samples,
            "Year": year,
            "Platform": gpl,
            "Summary": summary,
            "Queries": ", ".join(query_mappings.get(gse, []))
        })

print(f"\nFiltered to {len(parsed_data)} IBD anti-TNF datasets.")

# Output to JSON
with open("Code/geo_raw_results.json", "w", encoding="utf-8") as f:
    json.dump(parsed_data, f, indent=4)

# Print safely to terminal
print("\n=== RELEVANT GSE RECORDS ===")
for item in parsed_data:
    # Clean text for terminal print
    clean_title = item['Title'].encode('ascii', 'replace').decode('ascii')
    clean_summary = item['Summary'][:300].encode('ascii', 'replace').decode('ascii')
    print(f"Accession: {item['GSE']}")
    print(f"Title: {clean_title}")
    print(f"Samples: {item['n_samples']}")
    print(f"Year: {item['Year']}")
    print(f"Platform (GPL): {item['Platform']}")
    print(f"Queries: {item['Queries']}")
    print(f"Summary: {clean_summary}...")
    print("-" * 50)
