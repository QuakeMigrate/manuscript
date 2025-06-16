# -*- coding: utf-8 -*-
"""
This script downloads waveform data, response info & manual pick files for the Askja
QM-manpick comparison example presented in the manuscript:

    QuakeMigrate **

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
params = {'access_token': ACCESS_TOKEN}

INPUTS_DIR = pathlib.Path("./figures/Askja_QM-manpick_comparison/generate_results/inputs")

r = requests.get('https://zenodo.org/api/records/15236744/files',
                 params=params)

for repo_file in r.json():
    # download dataless directly -- does not need unzipping
    if repo_file["key"] == "dataless.xml":
        out_dir = INPUTS_DIR
        r_link = requests.get(repo_file["links"]["content"], params=params)
        with open(out_dir / "dataless.xml", "wb") as fid:
            fid.write(r_link.content)
    # unzip NLLoc obs files (pick files)
    elif repo_file["key"] == "obs.zip":
        out_dir = INPUTS_DIR / "NLLOC"
        # download, extract & save picks
        r_link = requests.get(repo_file["links"]["content"], params=params)
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)
    # unzip yearly waveform data dirs and write to mSEED dir
    else:
        out_dir = INPUTS_DIR / "mSEED"
        # download, extract & save waveform data
        r_link = requests.get(repo_file["links"]["content"], params=params)
        z = zipfile.ZipFile(io.BytesIO(r_link.content))
        z.extractall(out_dir)

# end
