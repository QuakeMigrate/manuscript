# -*- coding: utf-8 -*-
"""
This script runs Trigger for the Askja Volcanotectonic (VT) &
Deep-Long-Period (DLP) event example presented in the manuscript:

    QuakeMigrate **

"""

# Stop numpy using all available threads (these environment variables must be
# set before numpy is imported for the first time).
import os
os.environ.update(OMP_NUM_THREADS="1",
                  OPENBLAS_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1",
                  MKL_NUM_THREADS="1")

from quakemigrate import Trigger
from quakemigrate.io import read_lut

# --- i/o paths ---
lut_file = "./outputs/lut/askja.LUT"
run_path = "./outputs/runs"
run_name = "24h_run"

# --- Set time period over which to run trigger ---
starttime = "2011-10-26T00:00:00.0"
endtime = "2011-10-27T00:00:00.0"

# --- Load the LUT ---
lut = read_lut(lut_file=lut_file)

# --- Create new Trigger ---
trig = Trigger(lut, run_path=run_path, run_name=run_name, log=True,
               loglevel="info")

# --- Set trigger parameters ---
trig.marginal_window = 1.0
trig.min_event_interval = 2.0
trig.normalise_coalescence = True

# --- Static threshold ---
trig.threshold_method = "static"
trig.static_threshold = 1.27

# --- Toggle plotting options ---
trig.plot_trigger_summary = True
trig.xy_files = "./inputs/XY_FILES/askja_xyfiles.csv"

# --- Run trigger ---
trig.trigger(starttime, endtime, interactive_plot=False,
             region=[-17.1, 64.95, -3.0, -16.0, 65.30, 30.0])
