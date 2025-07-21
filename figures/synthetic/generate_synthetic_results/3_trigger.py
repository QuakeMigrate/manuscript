"""
This script runs the trigger stage for the synthetic example presented in the
manuscript:

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

from quakemigrate import Trigger
from quakemigrate.io import read_lut


# --- i/o paths ---
lut_file = "./outputs/lut/example.LUT"
run_path = "./outputs/runs"
run_name = "example_run"

# --- Set time period over which to run trigger ---
starttime = "2021-02-18T12:03:50.0"
endtime = "2021-02-18T12:06:10.0"

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)

# --- Create new Trigger ---
trig = Trigger(lut, run_path=run_path, run_name=run_name, log=True, loglevel="info")

# --- Set trigger parameters ---
trig.marginal_window = 0.2
trig.min_event_interval = 6.0
trig.normalise_coalescence = True

# --- Static threshold ---
trig.threshold_method = "static"
trig.static_threshold = 4.0

# --- Run trigger ---
trig.trigger(starttime, endtime, interactive_plot=False)
