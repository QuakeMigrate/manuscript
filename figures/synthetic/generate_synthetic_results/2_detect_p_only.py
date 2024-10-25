"""
This script runs the detect stage (P only) for the synthetic example presented in the
manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

from quakemigrate import QuakeScan
from quakemigrate.io import Archive, read_lut, read_stations
from quakemigrate.signal.onsets import STALTAOnset


# --- i/o paths ---
station_file = "./inputs/synthetic_stations.txt"
data_in = "./inputs/mSEED"
lut_out = "./outputs/lut/example.LUT"
run_path = "./outputs/runs"
run_name = "example_run_p_only"

# --- Set time period over which to run detect ---
starttime = "2021-02-18T12:03:50.0"
endtime = "2021-02-18T12:06:10.0"

# --- Read in station file ---
stations = read_stations(station_file)

# --- Create new Archive and set path structure ---
archive = Archive(
    archive_path=data_in, stations=stations, archive_format="YEAR/JD/STATION"
)

# --- Load the LUT ---
lut = read_lut(lut_file=lut_out)
lut.decimate([2, 2, 2], inplace=True)

# --- Create new Onset ---
onset = STALTAOnset(position="classic", sampling_rate=100)
onset.phases = ["P"]
onset.bandpass_filters = {"P": [1, 14, 2]}
onset.sta_lta_windows = {"P": [0.2, 1.5]}

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
scan.timestep = 120
scan.threads = 4  # NOTE: increase as your system allows to increase speed!

# --- Run detect ---
scan.detect(starttime, endtime)

