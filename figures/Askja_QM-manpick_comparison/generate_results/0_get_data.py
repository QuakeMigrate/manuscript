"""
This script downloads waveform data, response info & manual pick files for the Askja
QM-manpick comparison example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. Seismica, 2026.

    DOI: 10.26443/seismica.v5i1.1854

"""

import io
import pathlib
import re
import requests
import zipfile


INPUTS_DIR = pathlib.Path("./inputs")

MSEED_PATTERN = re.compile(r"^20...zip$")

r = requests.get("https://zenodo.org/api/records/15236744/files")

for repo_file in r.json()["entries"]:
    # download dataless directly -- does not need unzipping
    if repo_file["key"] == "dataless.xml":
        out_dir = INPUTS_DIR
        out_dir.mkdir(parents=True, exist_ok=True)
        r_link = requests.get(repo_file["links"]["content"])
        with open(out_dir / "dataless.xml", "wb") as fid:
            fid.write(r_link.content)
    # unzip NLLoc obs files (pick files)
    elif repo_file["key"] == "obs.zip":
        out_dir = INPUTS_DIR / "NLLOC"
        out_dir.mkdir(parents=True, exist_ok=True)
        # download, extract & save picks
        r_link = requests.get(repo_file["links"]["content"])
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)
    # unzip yearly waveform data dirs and write to mSEED dir
    elif MSEED_PATTERN.match(repo_file["key"]):
        out_dir = INPUTS_DIR / "mSEED"
        out_dir.mkdir(parents=True, exist_ok=True)
        # download, extract & save waveform data
        r_link = requests.get(repo_file["links"]["content"])
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)
