import urllib.request
import json
import time
import urllib.parse

queries = [
    "anti-TNF response Crohn's disease microarray",
    "infliximab response inflammatory bowel disease baseline",
    "adalimumab ulcerative colitis RNA-seq response",
    "anti-TNF non-response IBD gene expression"
]

def search_geo(query):
    print(f"Searching: {query}")
    encoded_query = urllib.parse.quote_plus(query)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={encoded_query}&retmode=json&retmax=50"
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
    id_str = ",".join(ids)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_str}&retmode=json"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            return data.get("result", {})
    except Exception as e:
        print(f"Error fetching summaries: {e}")
        return {}

all_results = {}
for q in queries:
    ids = search_geo(q)
    if ids:
        summaries = get_summary(ids)
        for uid in ids:
            if uid in summaries:
                all_results[uid] = summaries[uid]
        time.sleep(1) # Be nice to NCBI

print(f"\nFetched {len(all_results)} total unique GDS summaries.")

# Parse and print summaries
print(f"Sample info keys: {list(all_results.values())[0].keys() if all_results else 'No results'}")
if all_results:
    sample_key = list(all_results.keys())[0]
    if sample_key != 'uids':
        print(f"Sample info for {sample_key}: {list(all_results.values())[0]}")
    else:
        print(f"Sample info: {all_results}")

# Parse and print summaries
parsed_data = []
for uid, info in all_results.items():
    if uid == 'uids':
        continue
    entrytype = info.get("entrytype", "").upper()
    if entrytype != "GSE":
        continue
    gse = info.get("accession") or info.get("extid")
    title = info.get("title")
    summary = info.get("summary")
    n_samples = info.get("n_samples")
    pdate = info.get("pdat")
    year = pdate.split("/")[0] if pdate else "Unknown"
    gpl = info.get("gpl")
    
    parsed_data.append({
        "GSE": gse,
        "Title": title,
        "n_samples": n_samples,
        "Year": year,
        "Summary": summary,
        "GPL": gpl
    })

print(f"\nFound {len(parsed_data)} unique GSE records.")
print("\n=== GSE RECORDS ===")
for item in parsed_data:
    print(f"Accession: {item['GSE']}")
    print(f"Title: {item['Title']}")
    print(f"Samples: {item['n_samples']}")
    print(f"Year: {item['Year']}")
    print(f"Platform (GPL): {item['GPL']}")
    print(f"Summary: {item['Summary'][:300]}...")
    print("-" * 50)
