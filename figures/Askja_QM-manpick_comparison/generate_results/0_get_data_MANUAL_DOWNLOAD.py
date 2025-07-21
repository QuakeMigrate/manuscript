# -*- coding: utf-8 -*-
"""
This script is to be used where the user has manually downloaded the relevant
files from the Zenodo repository to a local directory -- it extracts waveform data,
response info & manual pick files and places them in the correct directories to run the
Askja QM-manpick comparison example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import zipfile
import pathlib


# Please ensure all files in the Zenodo repository (found at:
# https://zenodo.org/records/15236744 ) have been downloaded to a local directory. You
# will require a sharing link to gain access before the repository has been published.

# NOTE: set the "DOWNLOAD_DIR" variable to the location where the zip files and dataless
# inventory were saved. E.g. "/home/user/Downloads"
DOWNLOAD_DIR = pathlib.Path("/PATH/TO/DOWNLOAD_DIR")
INPUTS_DIR = pathlib.Path("./figures/Askja_QM-manpick_comparison/generate_results/inputs")


# Move dataless file to correct location
dataless_file = "dataless.xml"
(DOWNLOAD_DIR / dataless_file).rename(INPUTS_DIR / dataless_file)

# Extract NLLoc obs files (pick files)
out_dir = INPUTS_DIR / "NLLOC"
z = zipfile.ZipFile(DOWNLOAD_DIR / "obs.zip")
z.extractall(out_dir)

# Loop through waveform data directories & unzip
out_dir = INPUTS_DIR / "mSEED"
for year_zip in [f"{y}.zip" for y in range(2007, 2016)]:
    print(f"Working on {year_zip}...")
    z = zipfile.ZipFile(DOWNLOAD_DIR / year_zip)
    z.extractall(out_dir)

# end
