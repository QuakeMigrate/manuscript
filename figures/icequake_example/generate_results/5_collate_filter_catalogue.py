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
    pd.read_csv(event_file) for event_file in event_file_path.glob("*.event")
]

catalogue = pd.concat(all_event_dfs)


def product(row):
    return row["COV_ErrX"] * row["COV_ErrY"] * row["COV_ErrZ"]


def geometric_mean(row):
    return row["GlobalCovarianceProduct"] ** (1 / 3)


catalogue["GlobalCovarianceProduct"] = catalogue.apply(product, axis=1)
catalogue["GlobalCovarianceGeometricMean"] = catalogue.apply(geometric_mean, axis=1)

catalogue.to_csv("rutford_icequakes_all.csv", index=False)

# Filter
catalogue = catalogue[catalogue["COA"] >= 5.0]
catalogue = catalogue[catalogue["GlobalCovarianceGeometricMean"] <= 0.5]

catalogue.to_csv("rutford_icequakes_gc500_coa6.csv", index=False)
