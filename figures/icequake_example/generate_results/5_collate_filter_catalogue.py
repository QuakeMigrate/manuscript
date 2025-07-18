"""
This script collates all located events for the Rutford cryoseismicity example and
filters them, as set out in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import pandas as pd


event_file_path = pathlib.Path.cwd() / "outputs/runs/paper_run/locate/events"
all_event_dfs = [
    pd.read_csv(event_file) for event_file in sorted(event_file_path.glob("*.event"))
]

catalogue = pd.concat(all_event_dfs)

# Filter
filt_catalogue = catalogue.query("COV_Err_XYZ <= 0.15 and COA >= 5.5")

# Save for plotting
catalogue.to_csv("rutford_icequakes.csv", index=False)
filt_catalogue.to_csv("rutford_icequakes_gc150_coa55.csv", index=False)
