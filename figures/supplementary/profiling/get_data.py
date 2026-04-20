"""
This script downloads memray profile data for the computational cost profiling of the
Askja VT-DLP example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import io
import pathlib
import requests
import zipfile


PROFILES_DIR = pathlib.Path.cwd()

r = requests.get("https://zenodo.org/api/records/15236744/files")

for repo_file in r.json()["entries"]:
    # unzip profiles.zip
    if repo_file["key"] == "profiles.zip":
        out_dir = PROFILES_DIR / "profiles"
        out_dir.mkdir(parents=True, exist_ok=True)
        # download, extract & save picks
        r_link = requests.get(repo_file["links"]["content"])
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)
