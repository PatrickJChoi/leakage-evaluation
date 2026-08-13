import urllib.request
import re

url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
print(f"Listing platforms directories to find GPL570...")
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as response:
        html = response.read().decode('utf-8')
    hrefs = re.findall(r'href="([^"]+)"', html)
    unique_hrefs = sorted(list(set(hrefs)))
    for h in unique_hrefs:
        # If it doesn't contain a number >= 10, let's see. Let's print anything that starts with GPL and is short
        if h.startswith("GPL"):
            # Check if it matches GPL[0-9]nnn/
            match = re.match(r'GPL(\d+)?nnn/?', h)
            if match:
                # print matches
                print(f"Directory: {h}")
except Exception as e:
    print(f"Error: {e}")
