import urllib.request
import re

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as response:
        html = response.read().decode('utf-8')
    links = re.findall(r'href="([^"]+)"', html)
    print("Files in GSE120575 suppl directory:")
    for link in links:
        if not link.startswith('?') and not link.startswith('/') and not link.startswith('.'):
            print(" ", link)
except Exception as e:
    print(f"Error: {e}")
