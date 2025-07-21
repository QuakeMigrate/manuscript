"""
This script runs the trigger stage for the Askja volcano (Iceland)
Volcanotectonic (VT) & Deep-Long-Period (DLP) event example presented in the
manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

from quakemigrate import Trigger
from quakemigrate.io import read_lut

# --- i/o paths ---
lut_file = "./outputs/lut/askja.LUT"
run_path = "./outputs/runs"
run_name = "profiling"

# --- Set time period over which to run trigger ---
starttime = "2011-299T12:00:00.0"
endtime = "2011-299T13:00:00.0"

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)

# --- Create new Trigger ---
trig = Trigger(lut, run_path=run_path, run_name=run_name, log=True, loglevel="info")

# --- Set trigger parameters ---
trig.marginal_window = 1.0
trig.min_event_interval = 2.0
trig.normalise_coalescence = True

# --- Static threshold ---
trig.threshold_method = "static"
trig.static_threshold = 1.45

# --- Apply smoothing to coalescence trace before triggering ---
trig.smooth_coa = True
trig.smoothing_kernel_sigma = 0.25
trig.smoothing_kernel_width = 2

# --- Toggle plotting options ---
trig.plot_trigger_summary = True
trig.xy_files = "./inputs/XY_FILES/askja_xyfiles.csv"

# --- Run trigger ---
trig.trigger(
    starttime,
    endtime,
    interactive_plot=False,
    region=[-17.1, 64.95, -3.0, -16.0, 65.30, 30.0],
)
