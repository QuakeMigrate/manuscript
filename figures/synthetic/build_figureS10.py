"""
This script builds Figure S10 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from nllgrid import NLLGrid
from quakemigrate.io import read_lut


plt.style.use("../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"

lut = read_lut(lut_file="./generate_synthetic_results/outputs/lut/example.LUT")
grid = NLLGrid(f"generate_synthetic_results/time/layer.P.EQ.angle.buf")
dip = grid.dip[0]
earthquake_coords = [0.0, 0.0, 15.0]

# Extract grid axes
nx, nz = dip.shape
x0, y0, z0 = grid.x_orig, grid.y_orig, grid.z_orig
dx, dy, dz = grid.dx, grid.dy, grid.dz

# In GRID2D: x = range, y = depth
ranges = x0 + np.arange(nx) * dx
depths = z0 + np.arange(nz) * dz

# Depth axis should be positive downward
if depths[1] < depths[0]:
    depths = depths[::-1]
    dip = dip[:, ::-1]

fig, ax = plt.subplots(figsize=(7.08661, 7.08661 / 2), constrained_layout=True)

im = ax.pcolormesh(ranges, depths, dip.T, shading="auto", cmap="viridis")

ax.set_xlabel("Radial distance / km")
ax.set_ylabel("Depth / km")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"Dip angle / $\degree$")

ax.scatter(
    earthquake_coords[0],
    earthquake_coords[2],
    marker="*",
    zorder=12,
    s=50,
    c="#c51b8a",
    edgecolors="k",
    lw=0.35,
)

levels = np.arange(0, 190, 10)
cs = ax.contour(ranges, depths, dip.T, levels=levels, colors="white", linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=8, fmt="%d°")

earthquake_xyz = lut.coord2grid(earthquake_coords)[0]
for sx, sy, sz in lut.stations_xyz:
    dx, dy = sx - earthquake_xyz[0], sy - earthquake_xyz[1]
    r_offset = np.sqrt(dx**2 + dy**2)

    ax.scatter([r_offset], [sz], s=50, marker="^", color="black", zorder=10)

ax.set_ylim([30.0, -2.0])
ax.set_xlim([0.0, 25.0])

plt.savefig("figureS10.png", dpi=400)
