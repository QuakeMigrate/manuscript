"""
This script runs the trigger stage for the Rutford cryoseismicity example presented in
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

from quakemigrate import Trigger
from quakemigrate.io import read_lut


# --- i/o paths ---
lut_file = "./outputs/lut/icequake.LUT"
run_path = "./outputs/runs"
run_name = "paper_run"

# --- Set time period over which to run trigger ---
starttime = "2009-01-20T00:00:00.0"
endtime = "2009-01-21T00:00:00.0"

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)

# --- Create new Trigger ---
trig = Trigger(
    lut,
    run_path=run_path,
    run_name=run_name,
    # trigger_name="smoothing_on",
    output_all_columns=True,
    log=True,
    loglevel="info",
)

# --- Set trigger parameters ---
# For a complete list of parameters and guidance on how to choose them, please
# see the manual and read the docs.
trig.marginal_window = 0.06
trig.min_event_interval = 0.12
trig.normalise_coalescence = True

# --- Smoothing ---
# trig.smooth_coa = True
# trig.smoothing_kernel_sigma = 0.02
# trig.smoothing_kernel_width = 4

# --- Static threshold ---
trig.threshold_method = "static"
trig.static_threshold = 3.0

# --- Dynamic (Median Absolute Deviation) threshold ---
# trig.threshold_method = "mad"
# trig.mad_window_length = 300.0
# trig.mad_multiplier = 8.0

# --- Run trigger ---
trig.trigger(starttime, endtime, interactive_plot=False)
