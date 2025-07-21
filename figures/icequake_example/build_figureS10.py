"""
This script builds Supplementary Figure S10 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


plt.style.use("../../qm_manuscript.mplstyle")
plt.rcParams.update({"font.family": "Helvetica"})


# Read in raw catalogue
catalogue = pd.read_csv(pathlib.Path.cwd() / "generate_results/rutford_icequakes.csv")

# Set up figure
fig = plt.figure(figsize=(18 / 2.54, 8 / 2.54), facecolor="w")
axs = fig.subplots(ncols=2)

# Histogram for COV_Err_XYZ
catalogue.COV_Err_XYZ.hist(ax=axs[0], bins=np.arange(0, 1.8, 0.025), color="royalblue")
axs[0].axvline(0.15, c="r", ls="--", label="Filter value = 0.150 km")

axs[0].set_xlabel("COV_Err_XYZ / km")
axs[0].set_ylabel("Frequency")

axs[0].legend(fontsize=8)

# Histogram for COA
catalogue.COA.hist(
    ax=axs[1], bins=np.arange(0, 80, 0.2), label="Raw catalogue", color="royalblue"
)
catalogue.query("COV_Err_XYZ < 0.15").COA.hist(
    ax=axs[1],
    bins=np.arange(0, 80, 0.2),
    label="COV_Err_XYZ < 0.15 km",
    color="limegreen",
)
axs[1].axvline(5.5, c="r", ls="--", label="Filter value = 5.5")

axs[1].set_xlim(0, 30)

axs[1].set_xlabel("COA")
axs[1].set_ylabel("Frequency")

axs[1].legend(fontsize=8)


# Add label
for ax, letter in zip(axs, ["a", "b"]):
    ax.add_artist(
        AnchoredText(
            letter,
            loc="upper left",
            prop={"size": 9, "weight": "bold"},
            frameon=False,
            borderpad=0.1,
        )
    )

fig.tight_layout()
fig.savefig("./figureS10.png", dpi=400, bbox_inches="tight")
