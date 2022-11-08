# -*- coding: utf-8 -*-
"""
This script creates the lookup table and synthetic waveforms for the toy
example used in the manuscript:

    QuakeMigrate **

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from quakemigrate.io import read_lut


plt.style.use("qm_manuscript")
mpl.rcParams["font.family"] = "Helvetica"

lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")

# Make LUT figure
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7.08661, 7.08661/2), constrained_layout=True)

ax = axes[0]
ax.scatter(
    lut.station_data.Longitude,
    lut.station_data.Latitude,
    marker="^",
    color="black",
    s=20
)
ax.scatter(0.0, 0.0, c="#c51b8a", marker="*", s=35)

ax.set_xlabel("Longitude")
ax.set_xlim([-0.15, 0.15])

ax.set_ylabel("Latitude")
ax.set_ylim([-0.15, 0.15])

ax = axes[1]
ax.plot(lut.velocity_model.Vp, lut.velocity_model.Depth, c="#1b9e77", label="Vp")
ax.plot(lut.velocity_model.Vs, lut.velocity_model.Depth, c="#7570b3", label="Vs")

ax.legend(fontsize=8, loc=3, frameon=False)

ax.set_xlabel("Velocity / km s$^{-1}$")
ax.set_xlim([0, 8])

ax.set_ylabel("Depth / km")
ax.set_ylim(30, -1)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(1))

spine_weight = 1.2
for ax in axes:
    plt.setp(ax.spines.values(), linewidth=spine_weight)

    ax.xaxis.set_tick_params(width=spine_weight)
    ax.yaxis.set_tick_params(width=spine_weight)

    ax.tick_params(length=3)
    ax.tick_params(which="minor", length=2.25, width=1)

plt.savefig("figure3.pdf", dpi=400)
