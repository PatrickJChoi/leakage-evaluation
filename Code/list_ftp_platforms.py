import urllib.request
import re

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
print(f"Listing files in: {url}")
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as response:
        html = response.read().decode('utf-8')
    hrefs = re.findall(r'href="([^"]+)"', html)
    print("Files/Directories found (first 30):")
    # Sort and take first 30
    unique_hrefs = sorted(list(set(hrefs)))
    for h in unique_hrefs[:30]:
        if not h.startswith('/') and not h.startswith('?'):
            print(f"  - {h}")
except Exception as e:
    print(f"Error: {e}")
