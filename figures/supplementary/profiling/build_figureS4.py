"""
This script builds Figure S4 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from memray import FileReader


plt.style.use("qm_manuscript")
mpl.rcParams["font.family"] = "Helvetica"


def get_runtime(reader: FileReader) -> float:
    """Computes the runtime for a profile run."""

    return (reader.metadata.end_time - reader.metadata.start_time).total_seconds()


def get_high_watermark_in_gib(reader: FileReader) -> float:
    """Identify the peak memory usage and report it in GiB."""

    return reader.metadata.peak_memory / 1e9


profiles = {}
for profile in (pathlib.Path.cwd() / "profiles").glob("*.bin"):
    if "askja" in str(profile):
        continue
    reader = FileReader(profile)

    _, decimation_factor, n_threads, timestep = profile.stem.split(".")

    profiles[profile.stem] = {
        "reader": reader,
        "decimation_factor": decimation_factor,
        "n_threads": n_threads,
        "timestep": timestep,
    }


fig, axes = plt.subplots(
    ncols=3, nrows=2, figsize=(18.5 / 2.54, 8 / 2.54), constrained_layout=True
)

# Decimation factor
decimation_factors = [1, 2, 4]
runtimes, high_watermarks, average_memory = [], [], []
for decimation_factor in decimation_factors:
    profile = profiles[
        f"qm-detect.decimation-factor_{decimation_factor}.n-threads_16.timestep_300"
    ]
    runtimes.append(get_runtime(profile["reader"]))
    high_watermarks.append(get_high_watermark_in_gib(profile["reader"]))
    average_memory.append(
        np.mean(profile["reader"].get_temporal_high_water_mark_allocation_records()[1])
        / 1e9
    )

# Runtime
ax = axes[0][0]
ax.scatter(decimation_factors, runtimes, marker="s", c="#50aba3", s=10)
ax.set_xticks(decimation_factors)
ax.set_xticklabels([])
ax.set_ylabel("Runtime, s")
ax.set_ylim([0, 300])

# Decimation factor, memory usage
ax = axes[1][0]
ax.scatter(decimation_factors, high_watermarks, marker="d", c="#6a2764", s=10)
ax.scatter(decimation_factors, average_memory, marker="d", c="#a1d8e2", s=10)
ax.set_xticks(decimation_factors)
ax.set_xlabel("Decimation factor")
ax.set_ylabel("Memory high watermark / GiB")
ax.set_ylim([0, 25])

# Timestep
timesteps = [50, 300, 600, 1200, 1800]
runtimes, high_watermarks, average_memory = [], [], []
for timestep in timesteps:
    profile = profiles[
        f"qm-detect.decimation-factor_2.n-threads_16.timestep_{timestep}"
    ]
    runtimes.append(get_runtime(profile["reader"]))
    high_watermarks.append(get_high_watermark_in_gib(profile["reader"]))
    average_memory.append(
        np.mean(profile["reader"].get_temporal_high_water_mark_allocation_records()[1])
        / 1e9
    )

# Runtime
ax = axes[0][1]
ax.scatter(timesteps, runtimes, marker="s", c="#50aba3", s=10)
ax.set_xticks(timesteps)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_ylim([0, 300])

# Decimation factor, memory usage
ax = axes[1][1]
ax.scatter(timesteps, high_watermarks, marker="d", c="#6a2764", s=10)
ax.scatter(timesteps, average_memory, marker="d", c="#a1d8e2", s=10)
ax.set_xticks(timesteps)
ax.set_xlabel("Timestep / s")
ax.set_yticklabels([])
ax.set_ylim([0, 25])

# n_threads
n_threads = [1, 2, 4, 8, 16, 24, 32]
runtimes, high_watermarks, average_memory = [], [], []
for n_thread in n_threads:
    profile = profiles[
        f"qm-detect.decimation-factor_2.n-threads_{n_thread}.timestep_300"
    ]
    runtimes.append(get_runtime(profile["reader"]))
    high_watermarks.append(get_high_watermark_in_gib(profile["reader"]))
    average_memory.append(
        np.mean(profile["reader"].get_temporal_high_water_mark_allocation_records()[1])
        / 1e9
    )

# Runtime
ax = axes[0][2]
ax.scatter(n_threads, runtimes, marker="s", c="#50aba3", s=10)
ax.set_xticks(n_threads)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_ylim([0, 300])

# Decimation factor, memory usage
ax = axes[1][2]
ax.scatter(n_threads, high_watermarks, marker="d", c="#6a2764", s=10)
ax.scatter(n_threads, average_memory, marker="d", c="#a1d8e2", s=10)
ax.set_xticks(n_threads)
ax.set_xlabel("Number of threads")
ax.set_yticklabels([])
ax.set_ylim([0, 25])

fig.savefig("FigureS4.png", dpi=400)
