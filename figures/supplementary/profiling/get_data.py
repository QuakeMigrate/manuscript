"""
This script downloads memray profile data for the computational cost profiling of the
Askja VT-DLP example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import requests
import zipfile
import io
import os
import pathlib


# NOTE: to use the Zenodo API (reqiured for this script) you need to generate an access
# token for your Zenodo account; see https://zenodo.org/account/settings/applications/ .
# If this is not possible, you may also download directly through the web interface, and
# use the alternative "MANUAL_DOWNLOAD" script to extract the contents to the correct
# location. This will be the only available method until the Zenodo repository has been
# published.

# set your access token (from Zenodo) in the shell, with e.g. in Linux/bash shell:
# 'export ACCESS_TOKEN=12345'. If this is not possible, you can also manually set it
# below -- but do not accidentally upload to GitHub!
ACCESS_TOKEN = os.environ["ACCESS_TOKEN"]
params = {"access_token": ACCESS_TOKEN}

PROFILES_DIR = pathlib.Path("./figures/supplementary/profiling")

r = requests.get("https://zenodo.org/api/records/15236744/files", params=params)

for repo_file in r.json():
    # unzip profiles.zip
    if repo_file["key"] == "profiles.zip":
        out_dir = PROFILES_DIR / "profiles"
        # download, extract & save picks
        r_link = requests.get(repo_file["links"]["content"], params=params)
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)

# end
