import urllib.request
import tarfile
import os
import time

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE191nnn/GSE191328/suppl/GSE191328_RAW.tar"
dest = "GSE191328_RAW.tar"

baseline_gsms = [
    'GSM5743707', 'GSM5743715', 'GSM5743718', 'GSM5743723', 'GSM5743724', 
    'GSM5743729', 'GSM5743730', 'GSM5743735', 'GSM5743741', 'GSM5743745', 
    'GSM5743754', 'GSM5743760', 'GSM5743765', 'GSM5743841', 'GSM5743846', 
    'GSM5743848', 'GSM5743851', 'GSM5743856', 'GSM5743857', 'GSM5743858', 
    'GSM5743859', 'GSM5743860', 'GSM5743867', 'GSM5743869', 'GSM5743882', 
    'GSM5743885', 'GSM5743887', 'GSM5743889', 'GSM5743890', 'GSM5743891', 
    'GSM5743892', 'GSM5743893', 'GSM5743894', 'GSM5743896', 'GSM5743899'
]

def download_file(url, dest):
    print(f"Downloading {url} to {dest}...")
    start_time = time.time()
    try:
        # Custom request to NCBI with headers
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=180) as response, open(dest, 'wb') as out_file:
            # We can download in chunks and print progress
            total_size = int(response.getheader('Content-Length', 0))
            block_size = 1024 * 1024 * 5 # 5MB chunks
            downloaded = 0
            
            while True:
                buffer = response.read(block_size)
                if not buffer:
                    break
                downloaded += len(buffer)
                out_file.write(buffer)
                if total_size:
                    percent = (downloaded / total_size) * 100
                    print(f"Downloaded {downloaded / (1024*1024):.1f} MB / {total_size / (1024*1024):.1f} MB ({percent:.1f}%)")
            
        duration = time.time() - start_time
        print(f"Download completed in {duration:.1f} seconds.")
        return True
    except Exception as e:
        print(f"Error downloading: {e}")
        return False

def extract_specific_samples():
    print(f"Extracting {len(baseline_gsms)} samples from {dest}...")
    extract_dir = "Data/GSE191328_extracted"
    os.makedirs(extract_dir, exist_ok=True)
    
    try:
        with tarfile.open(dest, "r") as tar:
            members = tar.getmembers()
            # Filter members that match our baseline GSMs
            target_members = []
            for m in members:
                for gsm in baseline_gsms:
                    if gsm in m.name:
                        target_members.append(m)
                        break
            
            print(f"Found {len(target_members)} target files in the tar archive.")
            for idx, m in enumerate(target_members, 1):
                # Extract file
                tar.extract(m, path=extract_dir)
                if idx % 10 == 0 or idx == len(target_members):
                    print(f"Extracted {idx}/{len(target_members)} files...")
                    
        print(f"Extraction complete! Files saved to {extract_dir}")
        # Delete the large tar file to save disk space
        print(f"Deleting raw tar file {dest}...")
        os.remove(dest)
        print("Tar file deleted successfully.")
        return True
    except Exception as e:
        print(f"Error extracting tar: {e}")
        return False

if __name__ == "__main__":
    if not os.path.exists(dest) and not os.path.exists("Data/GSE191328_extracted"):
        download_file(url, dest)
    if os.path.exists(dest):
        extract_specific_samples()
