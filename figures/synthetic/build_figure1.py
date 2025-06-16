"""
This script builds Figure 1 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd
from quakemigrate.core import centred_sta_lta
from quakemigrate.io import read_lut


plt.style.use("qm_manuscript")
mpl.rcParams["font.family"] = "Helvetica"

run_path = pathlib.Path.cwd() / "generate_synthetic_results/outputs/runs/example_run"
lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")
decimation_factor = 4
lut.decimate([decimation_factor, decimation_factor, 1], inplace=True)
node_count = lut.node_count[0]
earthquake_coords = [0.0, 0.0, 15.0]
earthquake_ijk = lut.index2coord(earthquake_coords, inverse=True)

# Make figure and axes using semantic layout
fig = plt.figure(figsize=(18 / 2.54, 9.5 / 2.54), constrained_layout=True)
ax_dict = fig.subplot_mosaic("ABBC;DEEF", empty_sentinel="X")

ax = ax_dict["A"]
ax.axhspan(
    ymin=node_count / 2 - 1,
    ymax=node_count / 2,
    xmin=0,
    xmax=1 / node_count,
    color="#21908cff",
    lw=0.1,
)

ax = ax_dict["D"]
for xpos in range(node_count):
    for ypos in np.arange(0.5, node_count / 2 + 0.1, 1):
        ax.axhspan(
            ymin=ypos,
            ymax=ypos + 1,
            xmin=xpos / node_count,
            xmax=(xpos + 1) / node_count,
            color="#c1c1c1",
            lw=0.1,
        )

for xpos in np.arange(node_count / 2):
    ax.axhspan(
        ymin=-0.5,
        ymax=0.5,
        xmin=xpos / node_count,
        xmax=(xpos + 1) / node_count,
        color="#c1c1c1",
        lw=0.1,
    )

ax.axhspan(
    ymin=-0.5,
    ymax=0.5,
    xmin=(node_count - 1) / (2 * node_count),
    xmax=(node_count + 1) / (2 * node_count),
    color="#21908cff",
    lw=0.1,
)

for ax in [ax_dict["A"], ax_dict["D"]]:
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlim([-node_count / 2, node_count / 2])
    ax.set_ylim([-node_count / 2, node_count / 2])

    for pos in np.arange(-node_count / 2, node_count / 2 + 0.1, 1):
        ax.axvline(pos, ymin=-(node_count + 1), ymax=node_count, c="#e1e1e1", lw=0.5)
        ax.axhline(pos, xmin=-(node_count + 1), xmax=node_count, c="#e1e1e1", lw=0.5)

    ax.scatter(
        0.0, 0.0, marker="*", zorder=10, s=22, c="#c51b8a", edgecolors="k", lw=0.35
    )

    station_lons, station_lats, _ = lut.index2coord(
        list(
            zip(
                lut.station_data.Longitude,
                lut.station_data.Latitude,
                lut.station_data.Elevation,
            )
        ),
        inverse=True,
    ).T
    ax.scatter(
        station_lons - node_count // 2,
        station_lats - node_count // 2,
        marker="^",
        color="black",
        s=8,
        zorder=10,
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

bad_ttimes = pd.DataFrame(
    list(
        zip(
            lut.station_data["Name"].values,
            lut.traveltime_to("P", [0, 16, 32]),
            lut.traveltime_to("S", [0, 16, 32]),
        )
    )
)

good_ttimes = pd.DataFrame(
    list(
        zip(
            lut.station_data["Name"].values,
            lut.traveltime_to("P", earthquake_ijk),
            lut.traveltime_to("S", earthquake_ijk),
        )
    )
)

simulated_stream = obspy.read(
    "./generate_synthetic_results/inputs/mSEED/2021/049/STA*.m"
)

for ax, ttimes, o_time in zip(
    [ax_dict["B"], ax_dict["E"]], [bad_ttimes, good_ttimes], [29973, 30000]
):
    ticks = []
    i = 0
    for stat_id in station_names:
        stw = int(round(0.1 * 50))
        ltw = int(round(1.5 * 50))

        tr = simulated_stream.select(station=stat_id, component="Z")[0]
        onset = centred_sta_lta(tr.data**2, stw, ltw)

        ax.plot(onset / max(onset) + i, color="#F03B20")

        tr = simulated_stream.select(station=stat_id, component="[N,E]")[0]
        onset = centred_sta_lta(tr.data**2, stw, ltw)

        ax.plot(onset / max(onset) + i, color="#3182BD")

        P_ttime = ttimes[ttimes[0] == stat_id][1].values[0]
        S_ttime = ttimes[ttimes[0] == stat_id][2].values[0]

        ax.scatter(
            x=o_time + int(100 * P_ttime),
            y=i + 1.4,
            c="#F03B20",
            marker="v",
            s=10,
            edgecolors="k",
            lw=0.35,
        )
        ax.scatter(
            x=o_time + int(100 * S_ttime),
            y=i + 1.4,
            c="#3182BD",
            marker="v",
            s=10,
            edgecolors="k",
            lw=0.35,
        )

        ticks.append(i)
        i += 2.5

    ax.axvline(x=30000, ymin=0, ymax=i, c="#c1c1c1", ls="--", lw=1)
    ax.axvline(x=o_time, ymin=0, ymax=i, c="#49006a", ls="--", lw=1)

    ax.set_xlim([29850, 31250])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])

ax_dict["E"].set_xlabel(r"Time $\longrightarrow$")

ax = ax_dict["C"]

coa_data = np.load(run_path / "locate/coalescence_maps/20210218120500160.npy")
slice_ = coa_data[:, :, 32, 0]

bounds = np.stack(lut.get_grid_extent(cells=True), axis=-1)
gminx, gmaxx = bounds[0]
gminy, gmaxy = bounds[1]

nx, ny = [dim + 1 for dim in slice_.shape]
grid1, grid2 = np.mgrid[gminx : gmaxx : nx * 1j, gminy : gmaxy : ny * 1j]
sc = ax.pcolormesh(
    grid1,
    grid2,
    slice_,
    edgecolors="face",
    cmap="viridis",
    vmin=0,
    vmax=0.5 * coa_data.max(),
)

ax = ax_dict["F"]
slice_idxs = np.unravel_index(coa_data.argmax(), coa_data.shape)
slice_ = coa_data[:, :, slice_idxs[2], slice_idxs[3]]

nx, ny = [dim + 1 for dim in slice_.shape]
grid1, grid2 = np.mgrid[gminx : gmaxx : nx * 1j, gminy : gmaxy : ny * 1j]
sc = ax.pcolormesh(
    grid1,
    grid2,
    slice_,
    edgecolors="face",
    cmap="viridis",
    vmin=0,
    vmax=0.5 * coa_data.max(),
)

plt.savefig("figure1.png", dpi=400)
