"""
This script builds Figure S1 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from quakemigrate.io import read_lut


plt.style.use("../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"

lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")

# Make LUT figure
fig, axes = plt.subplots(
    nrows=1, ncols=2, figsize=(7.08661, 7.08661 / 2), constrained_layout=True
)

ax = axes[0]
ax.scatter(
    lut.station_data.Longitude,
    lut.station_data.Latitude,
    marker="^",
    color="black",
    s=20,
)
ax.scatter(0.0, 0.0, c="#c51b8a", marker="*", s=35)

ax.text(
    0.04,
    0.96,
    "a",
    ha="center",
    va="center",
    transform=ax.transAxes,
    fontweight="bold",
)

ax.set_xlabel(r"Longitude / $\degree$E")
ax.set_xlim([-0.15, 0.15])

ax.set_ylabel(r"Latitude / $\degree$N")
ax.set_ylim([-0.15, 0.15])

ax.tick_params(top=True, right=True)

ax = axes[1]
ax.plot(lut.velocity_model.Vp, lut.velocity_model.Depth, c="#F03B20", label="V$_p$")
ax.plot(lut.velocity_model.Vs, lut.velocity_model.Depth, c="#3182BD", label="V$_s$")

ax.text(
    0.04,
    0.96,
    "b",
    ha="center",
    va="center",
    transform=ax.transAxes,
    fontweight="bold",
)

ax.legend(fontsize=8, loc=3, frameon=False)

ax.set_xlabel("Velocity / km s$^{-1}$")
ax.set_xlim([0, 8])

ax.set_ylabel("Depth / km")
ax.set_ylim(30, -1)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

ax.tick_params(which="both", top=True, left=True)

ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(1))

plt.savefig("figureS1.png", dpi=400)
