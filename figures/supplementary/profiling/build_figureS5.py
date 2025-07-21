"""
This script builds Figure S5 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from memray import FileReader


plt.style.use("../../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"


def get_runtime(reader: FileReader) -> float:
    """Computes the runtime for a profile run."""

    return (reader.metadata.end_time - reader.metadata.start_time).total_seconds()


profile_reader = FileReader("profiles/askja-lut.bin")

runtime = get_runtime(profile_reader)
memory_usage = [
    v / 1e9 for v in profile_reader.get_temporal_high_water_mark_allocation_records()[1]
]

timestep = runtime / (len(memory_usage) - 1)
timestamps = np.arange(0, runtime + timestep, timestep)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.08661, 3), constrained_layout=True)

ax.scatter(timestamps, memory_usage, marker="+", c="k", s=20, zorder=2)
ax.plot(timestamps, memory_usage, lw=1, c="k", zorder=3)

ax.set_xlim([0, runtime])
ax.set_ylim(bottom=0)

ax.set_xlabel("Time / s")
ax.set_ylabel("Memory used / GiB")

fig.savefig("FigureS5.png", dpi=400)
