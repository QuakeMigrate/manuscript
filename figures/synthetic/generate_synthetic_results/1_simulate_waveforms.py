"""
This script creates the waveforms for the synthetic example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import numpy as np
from quakemigrate.io import read_lut

from simulate import GaussianDerivativeWavelet, simulate_waveforms


lut = read_lut("./outputs/lut/example.LUT")

mseed_output_dir = pathlib.Path.cwd() / "inputs/mSEED/2021/049"
mseed_output_dir.mkdir(parents=True, exist_ok=True)

# Calculate synthetic wavelets and migrate by calculated traveltimes
np.random.seed(4)  # Fix seed for reproducible results

# --- Build wavelet ---
frequency, sps, half_timespan = 4.0, 100, 300.0
wavelet = GaussianDerivativeWavelet(frequency, sps, half_timespan)

earthquake_coords = [0.0, 0.0, 15.0]
aoi = 80
magnitude = 2.2

simulated_stream = simulate_waveforms(
    wavelet, earthquake_coords, lut, magnitude=magnitude, angle_of_incidence=aoi
)

for tr in simulated_stream:
    fname = f"inputs/mSEED/2021/049/{tr.stats.station}_{tr.stats.component}.m"
    tr.write(fname, format="MSEED")
