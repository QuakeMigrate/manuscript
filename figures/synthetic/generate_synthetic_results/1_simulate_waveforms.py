# -*- coding: utf-8 -*-
"""
This script creates the lookup table and synthetic waveforms for the toy
example used in the manuscript:

    QuakeMigrate **

"""

import pathlib

import numpy as np
from quakemigrate.io import read_lut

import simulate as simulate


lut = read_lut("./outputs/lut/example.LUT")

mseed_output_dir = pathlib.Path.cwd() / "inputs/mSEED/2021/049"
mseed_output_dir.mkdir(parents=True, exist_ok=True)

# Calculate synthetic wavelets and migrate by calculated traveltimes
np.random.seed(4)  # Fix seed for reproducible results

# --- Build wavelet ---
frequency, sps, time_span = 4.0, 100, 300.0
wavelet = simulate.GaussianDerivativeWavelet(frequency, sps, time_span)

earthquake_coords = [0.0, 0.0, 15.0]
aoi = 80
magnitude = 2.2

simulated_stream = simulate.simulate_waveforms(
    wavelet, earthquake_coords, lut, magnitude=magnitude, angle_of_incidence=aoi
)

for tr in simulated_stream:
    fname = f"inputs/mSEED/2021/049/{tr.stats.station}_{tr.stats.component}.m"
    tr.write(fname, format="MSEED")
