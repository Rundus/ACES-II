import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import os

# URL containing the PNG files
# base_url = "http://tid.uio.no/plasma/aurora/skn4/5577/2022/20221120/ut17/"
base_url = "http://tid.uio.no/plasma/aurora/skn4/6300/2022/20221120/ut17/"

# Folder to save images
# save_dir = "/home/connor/Data/ROCKETS/ACESII/all_sky/skibotn/5577"
save_dir = "/home/connor/Data/ROCKETS/ACESII/all_sky/skibotn/6300"
os.makedirs(save_dir, exist_ok=True)

# Get webpage
response = requests.get(base_url)
response.raise_for_status()

# Parse HTML
soup = BeautifulSoup(response.text, "html.parser")

# Find all .png links
png_urls = []

for link in soup.find_all("a"):
    href = link.get("href")

    if href and href.lower().endswith(".png"):
        full_url = urljoin(base_url, href)
        png_urls.append(full_url)

print(f"Found {len(png_urls)} PNG files")

# Download each PNG
for url in png_urls:
    filename = os.path.basename(url)
    filepath = os.path.join(save_dir, filename)

    try:
        img_response = requests.get(url, stream=True)
        img_response.raise_for_status()

        with open(filepath, "wb") as f:
            for chunk in img_response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"Saved: {filename}")

    except Exception as e:
        print(f"Failed to download {filename}: {e}")

print("Done!")