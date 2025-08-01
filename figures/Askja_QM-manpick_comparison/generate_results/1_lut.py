# -*- coding: utf-8 -*-
"""
This script generates the traveltime look-up table (LUT) for the Askja volcano
(Iceland) Volcanotectonic (VT) & Deep-Long-Period (DLP) event example presented
in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

from obspy.core import AttribDict
from pyproj import Proj

from quakemigrate.io import read_stations, read_vmodel
from quakemigrate.lut import compute_traveltimes

# --- i/o paths ---
station_file = "./inputs/askja_QM-manpick_stations.txt"
vmodel_file = "./inputs/askja_vmodel.txt"
lut_out = "./outputs/lut/askja_QM-manpick_comparison.LUT"

# --- Read in the station information file ---
stations = read_stations(station_file)

# --- Read in the velocity model file ---
vmodel = read_vmodel(vmodel_file, comment="#")

# --- Define the input and grid projections ---
gproj = Proj(
    proj="lcc",
    units="km",
    lon_0=-16.6,
    lat_0=65.1,
    lat_1=64.9,
    lat_2=65.3,
    datum="WGS84",
    ellps="WGS84",
    no_defs=True,
)
cproj = Proj(proj="longlat", datum="WGS84", ellps="WGS84", no_defs=True)

# --- Define the grid specifications ---
# AttribDict behaves like a Python dict, but also has '.'-style access.
grid_spec = AttribDict()
grid_spec.ll_corner = [-17.3, 64.85, -3.0]
grid_spec.ur_corner = [-15.8, 65.4, 37.0]
grid_spec.node_spacing = [1.0, 1.0, 1.0]
grid_spec.grid_proj = gproj
grid_spec.coord_proj = cproj

# --- 1-D velocity model LUT generation (using NonLinLoc eikonal solver) ---
lut = compute_traveltimes(
    grid_spec,
    stations,
    method="1dnlloc",
    vmod=vmodel,
    phases=["P", "S"],
    log=True,
    save_file=lut_out,
)
