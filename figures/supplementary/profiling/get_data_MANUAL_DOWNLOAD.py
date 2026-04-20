"""
This script downloads memray profile data for the computational cost profiling of the
Askja VT-DLP example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. Seismica, 2026.

    DOI: 10.26443/seismica.v5i1.1854

"""

import zipfile
import pathlib


# Please ensure all files in the Zenodo repository (found at:
# https://zenodo.org/records/15236744 ) have been downloaded to a local directory.

# NOTE: set the "DOWNLOAD_DIR" variable to the location where the zip files and dataless
# inventory were saved. E.g. "/home/user/Downloads"
DOWNLOAD_DIR = pathlib.Path("/PATH/TO/DOWNLOAD_DIR")
PROFILES_DIR = pathlib.Path.cwd()


# Extract profiling data
out_dir = PROFILES_DIR / "profiles"
z = zipfile.ZipFile(DOWNLOAD_DIR / "profiles.zip")
z.extractall(out_dir)

# end
