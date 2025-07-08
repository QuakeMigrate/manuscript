"""
This script builds Figure S2 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import simulate as simulate


plt.style.use("qm_manuscript")
mpl.rcParams["font.family"] = "Helvetica"

# Calculate synthetic wavelets and migrate by calculated traveltimes
np.random.seed(10)  # Fix seed for reproducible results

# --- Build wavelet ---
frequency, sps, time_span = 4.0, 100, 300.0
wavelet = simulate.GaussianDerivativeWavelet(frequency, sps, time_span)

ttime_noise = np.random.normal(scale=0.05, size=1)
amp_noise = np.random.normal(scale=0.1, size=len(wavelet.data))

noisy_wavelet = np.roll(
    wavelet.data.copy() + amp_noise, int(wavelet.sps * ttime_noise[0])
)

fig, axes = plt.subplots(
    ncols=1, nrows=2, figsize=(7.08661, 4), constrained_layout=True
)

ax = axes[0]
ax.plot(wavelet.time, wavelet.data, c="k")
ax.axvline(0, c="#ac2335", ls="--")

ax.set_xticks([])

ax = axes[1]
ax.plot(wavelet.time, noisy_wavelet, c="k")
ax.axvline(0, c="#ac2335", ls="--", label="Simulated onset time")
ax.axvline(ttime_noise, c="#32af76", ls="--", label="Noise-adjusted onset time")
ax.legend(fontsize=8, loc=3, frameon=False)

ax.set_xlabel("Time / s")

labels = ["Gaussian-derivative wavelet", "Noise-adjusted wavelet"]
for ax, axlabel, label in zip(axes, "ab", labels):
    ax.set_xlim([-0.4, 0.8])
    ax.set_yticks([-1, 0, 1])

    ax.text(
        0.02,
        0.92,
        axlabel,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontweight="bold",
    )
    ax.text(0.99, 0.9, label, ha="right", va="center", transform=ax.transAxes)

fig.savefig("figureS2.png", dpi=400)
