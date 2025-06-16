"""
This script runs the locate stage for the Rutford cryoseismicity example presented in
the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
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

from quakemigrate import QuakeScan
from quakemigrate.io import Archive, read_lut, read_stations
from quakemigrate.signal.onsets import STALTAOnset
from quakemigrate.signal.pickers import GaussianPicker


# --- i/o paths ---
station_file = "./inputs/rutford_stations.txt"
data_in = "./inputs/mSEED"
lut_out = "./outputs/lut/icequake.LUT"
run_path = "./outputs/runs"
run_name = "paper_run"

# --- Set time period over which to run locate ---
starttime = "2009-01-20T00:00:00.0"
endtime = "2009-01-21T00:00:00.0"

# --- Read in station file ---
stations = read_stations(station_file)

# --- Create new Archive and set path structure ---
archive = Archive(
    archive_path=data_in, stations=stations, archive_format="YEAR/JD/STATION"
)

# --- Load the LUT ---
lut = read_lut(lut_file=lut_out)

# --- Create new Onset ---
onset = STALTAOnset(position="centred", sampling_rate=250)
onset.phases = ["P", "S"]
onset.bandpass_filters = {"P": [20, 124, 4], "S": [10, 124, 4]}
onset.sta_lta_windows = {"P": [0.01, 0.25], "S": [0.05, 0.5]}
onset.channel_maps = {"P": "*[1]", "S": "*[2,3]"}

# --- Create new PhasePicker ---
picker = GaussianPicker(onset=onset)
picker.plot_picks = False

# --- Create new QuakeScan ---
scan = QuakeScan(
    archive,
    lut,
    onset=onset,
    picker=picker,
    run_path=run_path,
    run_name=run_name,
    log=True,
    loglevel="info",
)

# --- Set locate parameters ---
# For a complete list of parameters and guidance on how to choose them, please
# see the manual and read the docs.
scan.marginal_window = 0.06
scan.threads = 4  # NOTE: increase as your system allows to increase speed!

# --- Toggle plotting options ---
scan.plot_event_summary = False

# --- Toggle writing of waveforms ---
scan.write_cut_waveforms = True

# --- Run locate ---
scan.locate(starttime=starttime, endtime=endtime)
