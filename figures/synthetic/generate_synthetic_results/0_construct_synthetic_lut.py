"""
This script creates the lookup table for the synthetic example presented in the
manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import numpy as np
from obspy.core import AttribDict
import pandas as pd
from pyproj import Proj
from quakemigrate.io import read_vmodel
from quakemigrate.lut import compute_traveltimes


# Build synthetic lookup table
station_file = "./inputs/synthetic_stations.txt"
vmodel_file = "./inputs/velocity_model.csv"
lut_out = "./outputs/lut/example.LUT"

# --- Build station file ---
rng = np.random.default_rng(13)  # Fix seed for reproducible results
stations = pd.DataFrame()
stations["Network"] = ["SC"] * 10
stations["Name"] = [f"STA{i}" for i in range(10)]
stations["Longitude"] = rng.uniform(low=-0.15, high=0.15, size=10)
stations["Latitude"] = rng.uniform(low=-0.15, high=0.15, size=10)
stations["Elevation"] = rng.uniform(low=-0.0, high=1.0, size=10)
stations.to_csv(station_file, index=False)

# --- Read in the velocity model file ---
vmodel = read_vmodel(vmodel_file)

# --- Define the input and grid projections ---
gproj = Proj(
    proj="lcc",
    units="km",
    lon_0=0.0,
    lat_0=0.0,
    lat_1=-0.10,
    lat_2=0.101,
    datum="WGS84",
    ellps="WGS84",
    no_defs=True,
)
cproj = Proj(proj="longlat", datum="WGS84", ellps="WGS84", no_defs=True)

# --- Define the grid specifications ---
# AttribDict behaves like a Python dict, but also has '.'-style access.
grid_spec = AttribDict()
grid_spec.ll_corner = [-0.15, -0.15, -1.0]
grid_spec.ur_corner = [0.15, 0.15, 30.0]
grid_spec.node_spacing = [0.5, 0.5, 0.5]
grid_spec.grid_proj = gproj
grid_spec.coord_proj = cproj

# --- Homogeneous LUT generation ---
lut = compute_traveltimes(
    grid_spec,
    stations,
    method="1dnlloc",
    vmod=vmodel,
    phases=["P", "S"],
    log=True,
    save_file=lut_out,
)
print()
print(lut)
