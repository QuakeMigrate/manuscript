"""
This script builds Figure 5 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd
from matplotlib.patches import Ellipse, Rectangle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from quakemigrate.io import read_lut, read_coalescence


plt.style.use("../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"

lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")

run_path = pathlib.Path.cwd() / "generate_synthetic_results/outputs/runs/example_run"

marginalised_coa_map = read_coalescence(
    run_path / "locate/marginalised_coalescence_maps/20210218120500160.npy"
)

event = pd.read_csv(run_path / "locate/events/20210218120500160.event")

# Extract indices and grid coordinates of maximum coalescence
coa_map = np.ma.masked_invalid(marginalised_coa_map)
idx_max = np.column_stack(np.where(coa_map == np.nanmax(coa_map)))[0]
slices = [
    coa_map[:, :, idx_max[2]],
    coa_map[:, idx_max[1], :],
    coa_map[idx_max[0], :, :].T,
]
otime = obspy.UTCDateTime(event["DT"].values[0])

fig = plt.figure(figsize=(3.6, 5.5))

# --- Plot LUT, waveform gather, and max coalescence trace ---
station_clr = "white"
hypocentre = [event["X"].values[0], event["Y"].values[0], event["Z"].values[0]]

xy = plt.subplot2grid((7, 5), (0, 0), colspan=5, rowspan=5, fig=fig)
xz = plt.subplot2grid((7, 5), (5, 0), colspan=5, rowspan=2, fig=fig)

xz.get_shared_x_axes().joined(xy, xz)

# --- Set aspect ratio ---
# Aspect is defined such that a circle will be stretched so that its
# height is aspect times the width.
cells_extent = lut.get_grid_extent(cells=True)
extent = abs(cells_extent[1] - cells_extent[0])
grid_size = lut.node_spacing * lut.node_count
aspect = (extent[0] * grid_size[1]) / (extent[1] * grid_size[0])
xy.set_aspect(aspect=aspect)

bounds = np.stack(cells_extent, axis=-1)
for i, j, ax in [(0, 1, xy), (0, 2, xz)]:
    gminx, gmaxx = bounds[i]
    gminy, gmaxy = bounds[j]

    ax.set_xlim([gminx, gmaxx])
    ax.set_ylim([gminy, gmaxy])

    # --- Plot crosshair for event hypocentre ---
    if hypocentre is not None:
        ax.axvline(x=hypocentre[i], ls="--", lw=0.8, c="white")
        ax.axhline(y=hypocentre[j], ls="--", lw=0.8, c="white")

    slice_ = slices[i + j - 1]
    nx, ny = [dim + 1 for dim in slice_.shape]
    grid1, grid2 = np.mgrid[gminx : gmaxx : nx * 1j, gminy : gmaxy : ny * 1j]
    sc = ax.pcolormesh(grid1, grid2, slice_, edgecolors="face", cmap="viridis")

# # # --- Add colourbar ---
# cbax = fig.add_axes([0.68, 0.05, 0.22, 0.02])
# cb = fig.colorbar(sc, cax=cbax, orientation="horizontal")
# # cb.ax.set_xlabel("Normalised coalescence\nvalue", rotation=0, fontsize=5)

# --- Add colourbar ---
xz.add_patch(
    Rectangle(
        (0.62, 0.0),
        0.38,
        0.4,
        facecolor="w",
        edgecolor=None,
        alpha=0.9,
        transform=xz.transAxes,
    )
)

cax = inset_axes(
    xz,
    width="30%",
    height="7%",
    bbox_to_anchor=(0, 0.1, 0.99, 1),
    bbox_transform=xz.transAxes,
    loc="lower right",
)
fig.colorbar(sc, cax=cax, orientation="horizontal")
cax.set_xlabel("Normalised coalescence\nvalue", rotation=0, fontsize=5)
cax.xaxis.set_label_position("top")

xy.scatter(
    lut.station_data.Longitude.values,
    lut.station_data.Latitude.values,
    s=12,
    marker="^",
    zorder=20,
    c="white",
)
xz.scatter(
    lut.station_data.Longitude.values,
    lut.station_data.Elevation.values,
    s=8,
    marker="^",
    zorder=20,
    c="white",
)

# hypocentre = [pos.values[0] for pos in [event["GAU_X"], event["GAU_Y"], event["GAU_Z"]]]
error = [
    event["GAU_ErrX"].values[0],
    event["GAU_ErrY"].values[0],
    event["GAU_ErrZ"].values[0],
]
xyz = lut.coord2grid(hypocentre)[0]
d = abs(hypocentre - lut.coord2grid(xyz + error, inverse=True))[0]

xy_err = Ellipse(
    (hypocentre[0], hypocentre[1]),
    2 * d[0],
    2 * d[1],
    lw=1,
    edgecolor="k",
    fill=False,
    label="Gaussian uncertainty",
)
xz_err = Ellipse(
    (hypocentre[0], hypocentre[2]), 2 * d[0], 2 * d[2], lw=1, edgecolor="k", fill=False
)
xy.add_patch(xy_err)
xz.add_patch(xz_err)

# --- Add scalebar ---
num_cells = np.ceil(lut.node_count[0] / 10)
length = num_cells * lut.node_spacing[0]
size = extent[0] * length / grid_size[0]
scalebar = AnchoredSizeBar(
    xy.transData,
    size=size,
    label=f"{length} {lut.unit_name}",
    loc="lower right",
    pad=0.5,
    sep=5,
    frameon=False,
    color="white",
    fontproperties=fm.FontProperties(size=6),
)
xy.add_artist(scalebar)

# --- Axes labelling ---
xy.tick_params(
    which="both",
    left=True,
    right=True,
    top=True,
    bottom=True,
    labelleft=True,
    labeltop=False,
    labelright=False,
    labelbottom=False,
)
xy.set_ylabel(r"Latitude / $\degree$N")
xy.yaxis.set_label_position("left")

xy.text(
    0.03,
    0.03,
    "a",
    ha="center",
    va="center",
    transform=xy.transAxes,
    fontweight="bold",
    color="white",
)
xz.text(
    0.03,
    0.07,
    "b",
    ha="center",
    va="center",
    transform=xz.transAxes,
    fontweight="bold",
    color="white",
)

xz.invert_yaxis()
xz.tick_params(
    which="both",
    left=True,
    right=True,
    top=True,
    bottom=True,
    labelleft=True,
    labeltop=False,
    labelright=False,
    labelbottom=True,
)
xz.set_xlabel(r"Longitude / $\degree$E")
xz.set_ylabel(f"Depth / {lut.unit_name}")
xz.yaxis.set_label_position("left")

plt.savefig("figure5.png", dpi=400, bbox_inches="tight")
