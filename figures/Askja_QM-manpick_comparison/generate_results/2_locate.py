# -*- coding: utf-8 -*-
"""
This script runs locate from a mock trigger file containing a list of
manually picked earthquakes from the region around Askja volcano (Iceland)
for the purpose of benchmarking the location performance of QuakeMigrate, as
presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

# Stop numpy using all available threads (these environment variables must be
# set before numpy is imported for the first time).
import os

os.environ.update(
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
)

from obspy.core import AttribDict

from quakemigrate import QuakeScan
from quakemigrate.io import Archive, read_lut, read_stations, read_response_inv
from quakemigrate.signal.onsets import STALTAOnset
from quakemigrate.signal.pickers import GaussianPicker
from quakemigrate.signal.local_mag import LocalMag

# --- i/o paths ---
station_file = "./inputs/askja_QM-manpick_stations.txt"
response_file = "./inputs/dataless.xml"
data_in = "./inputs/mSEED"
lut_file = "./outputs/lut/askja_QM-manpick_comparison.LUT"
run_path = "./outputs/runs"
run_name = "QM-manpick_comparison"
trigger_file = "./inputs/manpick_events.csv"

# --- Read in station file ---
stations = read_stations(station_file)

# --- Read in response inventory ---
response_inv = read_response_inv(response_file)

# --- Specify parameters for response removal ---
response_params = AttribDict()
response_params.pre_filt = (0.05, 0.06, 20, 23)
response_params.water_level = 60
response_params.remove_full_response = False

# --- Create new Archive and set path structure ---
archive = Archive(
    archive_path=data_in,
    stations=stations,
    response_inv=response_inv,
    response_removal_params=response_params,
    resample=True,
    upfactor=5,
)
archive.format = "{year}/{jday:03d}/*"

# --- Specify parameters for amplitude measurement ---
amp_params = AttribDict()
amp_params.signal_window = 1.0
amp_params.noise_window = 5.0
amp_params.noise_measure = "ENV"
amp_params.bandpass_filter = True
amp_params.bandpass_lowcut = 2.0
amp_params.bandpass_highcut = 20.0
amp_params.filter_corners = 4


# --- Specify parameters for magnitude calculation ---
mag_params = AttribDict()
mag_params.A0 = "Greenfield2018_askja"
mag_params.use_hyp_dist = True
mag_params.amp_feature = "S_amp"
mag_params.trace_filter = ".*H[NE]$"
mag_params.noise_filter = 3.0

mags = LocalMag(amp_params=amp_params, mag_params=mag_params, plot_amplitudes=True)

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)

# --- Create new Onset ---
onset = STALTAOnset(
    position="centred", sampling_rate=50, signal_transform="env_squared"
)
onset.phases = ["P", "S"]
onset.bandpass_filters = {"P": [2, 16, 2], "S": [2, 14, 2]}
onset.sta_lta_windows = {"P": [0.2, 1.0], "S": [0.2, 1.0]}
onset.all_channels = False

# --- Create new PhasePicker ---
picker = GaussianPicker(onset=onset)
picker.plot_picks = False

# --- Create new QuakeScan ---
scan = QuakeScan(
    archive,
    lut,
    onset=onset,
    picker=picker,
    mags=mags,
    run_path=run_path,
    run_name=run_name,
    log=True,
    loglevel="info",
)

# --- Set locate parameters ---
scan.marginal_window = 1.0
scan.threads = 4  # NOTE: increase as your system allows to increase speed!

# --- Toggle plotting options ---
scan.plot_event_summary = True
scan.plot_all_stns = False
scan.xy_files = "./inputs/XY_FILES/askja_xyfiles.csv"

# --- Toggle writing of waveforms ---
scan.write_cut_waveforms = False

# --- Toggle writing of marginal coalescence ---
scan.write_marginal_coalescence = False

# --- Run locate ---
scan.locate(trigger_file=trigger_file)
