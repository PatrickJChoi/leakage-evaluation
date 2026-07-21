import urllib.request
import re

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL570/"
print(f"Listing files in: {url}")
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as response:
        html = response.read().decode('utf-8')
    hrefs = re.findall(r'href="([^"]+)"', html)
    print("Files/Directories found:")
    for h in set(hrefs):
        if not h.startswith('/') and not h.startswith('?'):
            print(f"  - {h}")
except Exception as e:
    print(f"Error: {e}")
