"""
This script runs the detect stage for the Askja volcano (Iceland)
Volcanotectonic (VT) & Deep-Long-Period (DLP) event example presented in the
manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

from quakemigrate import QuakeScan
from quakemigrate.io import Archive, read_lut, read_stations
from quakemigrate.signal.onsets import STALTAOnset

# --- i/o paths ---
station_file = "./inputs/askja_stations.txt"
data_in = "./inputs/mSEED"
lut_file = "./outputs/lut/askja.LUT"
run_path = "./outputs/runs"
run_name = "paper_run"

# --- Set time period over which to run detect ---
starttime = "2011-299T00:00:00.0"
endtime = "2011-300T00:00:00.0"

# --- Read in station file ---
stations = read_stations(station_file)

# --- Create new Archive and set path structure ---
archive = Archive(
    archive_path=data_in, stations=stations, archive_format="YEAR/JD/STATION"
)

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)
lut.decimate([2, 2, 2], inplace=True)

# --- Create new Onset ---
onset = STALTAOnset(
    position="classic", sampling_rate=50, signal_transform="env_squared"
)
onset.phases = ["P", "S"]
onset.bandpass_filters = {"P": [2, 16, 2], "S": [2, 14, 2]}
onset.sta_lta_windows = {"P": [0.2, 1.0], "S": [0.2, 1.0]}

# --- Create new QuakeScan ---
scan = QuakeScan(
    archive,
    lut,
    onset=onset,
    run_path=run_path,
    run_name=run_name,
    log=True,
    loglevel="info",
)

# --- Set detect parameters ---
scan.timestep = 300.0
scan.threads = 32  # NOTE: increase as your system allows to increase speed!

# --- Run detect ---
scan.detect(starttime, endtime)
