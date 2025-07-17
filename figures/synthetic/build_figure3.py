"""
This script builds Figure 3 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import obspy
from quakemigrate.core import centred_sta_lta


plt.style.use("../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"


def make_colours():
    return iter(plt.cm.viridis(np.linspace(0, 10, 11) % 10 / 10))


fig, axes = plt.subplots(
    nrows=2, ncols=2, figsize=(7.08661, 6.8), constrained_layout=True
)

colours = make_colours()
station_names = [
    "STA7",
    "STA8",
    "STA5",
    "STA3",
    "STA6",
    "STA0",
    "STA9",
    "STA1",
    "STA2",
    "STA4",
]
ticks = []

simulated_stream = obspy.read(
    "./generate_synthetic_results/inputs/mSEED/2021/049/STA*.m"
)

stw = int(round(0.1 * 50))
ltw = int(round(1.5 * 50))
i = 0
for stat_id, clr in zip(station_names, colours):
    tr = simulated_stream.select(station=stat_id, component="Z")[0]
    onset = centred_sta_lta(tr.data**2, stw, ltw)

    axes[0][0].plot(tr.data + i, color=clr)
    axes[0][1].plot(onset / max(onset) + i, color=clr)

    tr = simulated_stream.select(station=stat_id, component="[N,E]")[0]
    onset = centred_sta_lta(tr.data**2, stw, ltw)

    axes[1][0].plot(tr.data + i, color=clr)
    axes[1][1].plot(onset / max(onset) + i, color=clr)

    ticks.append(i)
    i += 2.5

for ax in axes.flatten():
    ax.set_xlim([30225, 30875])
    ax.set_ylim([-2, len(simulated_stream) * 2.5 / 3])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])

labels = [
    "Vertical component (P-wave)",
    "P-wave Onset",
    "Horizontal component (S-wave)",
    "S-wave Onset",
]
for ax, panel_label, label in zip(axes.flatten(), "abcd", labels):
    ax.text(
        0.03,
        0.96,
        panel_label,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontweight="bold",
    )
    ax.text(0.985, 0.96, label, ha="right", va="center", transform=ax.transAxes)

for ax in [axes.flatten()[0], axes.flatten()[2]]:
    ax.set_yticks(ticks)
    ax.set_yticklabels(
        station_names, size=7, fontname="monospace", color="w", fontweight="bold"
    )

    colours = make_colours()
    # Colour y labels to match lines
    gytl = ax.get_yticklabels()
    for yt, color in zip(gytl, colours):
        yt.set_bbox(dict(facecolor=color, edgecolor=color, pad=1))

axes[1][1].set_xlabel("T", c="white")
fig.suptitle(r"Time $\longrightarrow$", x=0.525, y=0.025)

plt.savefig("figure3.png", dpi=400)
