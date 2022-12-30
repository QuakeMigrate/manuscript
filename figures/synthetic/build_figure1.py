# -*- coding: utf-8 -*-
"""
This script builds Figure 1 of the manuscript:

    QuakeMigrate **

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import obspy
import pandas as pd

from quakemigrate.io import read_lut

import generate_synthetic_results.simulate as simulate


plt.style.use("qm_manuscript")
mpl.rcParams["font.family"] = "Helvetica"

lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")
decimation_factor = 4
lut.decimate(
    [decimation_factor, decimation_factor, 1],
    inplace=True
)
node_count = lut.node_count[0]
earthquake_coords = [0.0, 0.0, 15.0]
earthquake_ijk = lut.index2coord(earthquake_coords, inverse=True)

# Make figure and axes using semantic layout
fig = plt.figure(figsize=(7.08661, 3.8), constrained_layout=True)
ax_dict = fig.subplot_mosaic("BBAC;EEDF", empty_sentinel="X")

ax = ax_dict["A"]
ax.axhspan(
    ymin=node_count/2 - 1,
    ymax=node_count/2,
    xmin=0,
    xmax=1/node_count,
    color="#dd3497",
    lw=0.1
)

ax = ax_dict["D"]
for xpos in range(node_count):
    for ypos in np.arange(0.5, node_count/2+0.1, 1):
        ax.axhspan(
            ymin=ypos,
            ymax=ypos+1,
            xmin=xpos/node_count,
            xmax=(xpos+1)/node_count,
            color="#dd3497",
            lw=0.1
        )

for xpos in np.arange(node_count/2):
    ax.axhspan(
        ymin=-0.5,
        ymax=0.5,
        xmin=xpos/node_count,
        xmax=(xpos+1)/node_count,
        color="#dd3497",
        lw=0.1
    )

for ax in [ax_dict["A"], ax_dict["D"]]:
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlim([
        -node_count/2,
        node_count/2
    ])
    ax.set_ylim([
        -node_count/2,
        node_count/2
    ])
    
    for pos in np.arange(-node_count/2, node_count/2 + 0.1, 1):
        ax.axvline(pos, ymin=-(node_count+1), ymax=node_count, c="#e1e1e1", lw=0.5)
        ax.axhline(pos, xmin=-(node_count+1), xmax=node_count, c="#e1e1e1", lw=0.5)

    ax.scatter(
        0.,
        0.,
        marker="*",
        zorder=10,
        s=22,
        c="#c51b8a",
        edgecolors="k",
        lw=0.35
    )

    station_lons, station_lats, _ = lut.index2coord(
        list(zip(
            lut.station_data.Longitude,
            lut.station_data.Latitude,
            lut.station_data.Elevation
        )),
        inverse=True
    ).T
    ax.scatter(
        station_lons-node_count//2,
        station_lats-node_count//2,
        marker="^",
        color="black",
        s=8,
        zorder=10
    )

for ax in [ax_dict["C"], ax_dict["F"]]:
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])

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

bad_ttimes = pd.DataFrame(list(zip(
    lut.station_data["Name"].values,
    lut.traveltime_to("P", [0, 16, 32]),
    lut.traveltime_to("S", [0, 16, 32])
)))

good_ttimes = pd.DataFrame(list(zip(
    lut.station_data["Name"].values,
    lut.traveltime_to("P", earthquake_ijk),
    lut.traveltime_to("S", earthquake_ijk)
)))

simulated_stream = obspy.read(
    "./generate_synthetic_results/inputs/mSEED/2021/049/STA*.m"
)

for ax, ttimes, o_time in zip(
    [ax_dict["B"], ax_dict["E"]],
    [bad_ttimes, good_ttimes],
    [29958, 30058]
):
    ticks = []
    i = 0
    for stat_id in station_names:
        tr = simulated_stream.select(station=stat_id, component="Z")[0]
        onset, _ = simulate.sta_lta_onset(tr.data, fs=tr.stats.sampling_rate)

        ax.plot(onset / max(onset) + i, color="#F03B20")

        tr = simulated_stream.select(station=stat_id, component="[N,E]")[0]
        onset, _ = simulate.sta_lta_onset(tr.data, fs=tr.stats.sampling_rate)

        ax.plot(onset / max(onset) + i, color="#3182BD")

        P_ttime = ttimes[ttimes[0] == stat_id][1].values[0]
        S_ttime = ttimes[ttimes[0] == stat_id][2].values[0]

        ax.scatter(
            x=o_time+int(100*P_ttime),
            y=i+1.4,
            c="#F03B20",
            marker="v",
            s=10,
            edgecolors="k",
            lw=0.35
        )
        ax.scatter(
            x=o_time+int(100*S_ttime),
            y=i+1.4,
            c="#3182BD",
            marker="v",
            s=10,
            edgecolors="k",
            lw=0.35
        )

        ticks.append(i)
        i += 2.5
    
    ax.axvline(x=30071, ymin=0, ymax=i, c="#c1c1c1", ls="--", lw=1)
    ax.axvline(x=o_time+13, ymin=0, ymax=i, c="#49006a", ls="--", lw=1)
    
    ax.set_xlim([29900, 31250])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])

ax_dict["E"].set_xlabel("Time $\longrightarrow$")

coa_data = np.load(
    "generate_synthetic_results/4D_coalescence/4d_coa_2021-02-18T12:03:38.900000Z"
)

ax = ax_dict["C"]
slice_ = coa_data[:, :, 32, 3475]
ax.imshow(
    slice_,
    interpolation="spline36",
    cmap="plasma",
    vmin=0,
    vmax=10
)

ax = ax_dict["F"]
slice_ = coa_data[:, :, 32, 3525]
ax.imshow(
    slice_,
    interpolation="spline36",
    cmap="plasma",
    vmin=0,
    vmax=111.459
)

plt.savefig("figure1.png", dpi=400)
