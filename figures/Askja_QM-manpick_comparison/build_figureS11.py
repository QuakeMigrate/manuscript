"""
This script builds Supplementary Figure S11 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


plt.style.use("../../qm_manuscript.mplstyle")
plt.rcParams.update({"font.family": "Helvetica"})

# Read in catalogue
man_qm_combined = pd.read_csv("./generate_results/outputs/QM-man_combined.csv")

# Calculate ratio between quoted uncertainty and location difference - NLLoc uncertainties
man_qm_combined.loc[:, "delta_sigma_ratio_X_MAN"] = (
    man_qm_combined.X_diff.abs() / man_qm_combined.Err_X
)
man_qm_combined.loc[:, "delta_sigma_ratio_Y_MAN"] = (
    man_qm_combined.Y_diff.abs() / man_qm_combined.Err_Y
)
man_qm_combined.loc[:, "delta_sigma_ratio_Z_MAN"] = (
    man_qm_combined.Z_diff.abs() / man_qm_combined.Err_Z
)
# Calculate ratio between quoted uncertainty and location difference - QM uncertainties
man_qm_combined.loc[:, "delta_sigma_ratio_X_QM"] = (
    man_qm_combined.X_diff.abs() / man_qm_combined.GAU_ErrX
)
man_qm_combined.loc[:, "delta_sigma_ratio_Y_QM"] = (
    man_qm_combined.Y_diff.abs() / man_qm_combined.GAU_ErrY
)
man_qm_combined.loc[:, "delta_sigma_ratio_Z_QM"] = (
    man_qm_combined.Z_diff.abs() / man_qm_combined.GAU_ErrZ
)

# Plot figure
fig = plt.figure(figsize=(18 / 2.54, 7 / 2.54), facecolor="w")
axs = fig.subplots(ncols=3, sharey=True)

for ax, ordinal in zip(axs, ["X", "Y", "Z"]):
    sns.ecdfplot(
        data=man_qm_combined,
        x=f"delta_sigma_ratio_{ordinal}_QM",
        ax=ax,
        c=plt.cm.viridis(0.4),
        label="QuakeMigrate",
    )
    sns.ecdfplot(
        data=man_qm_combined,
        x=f"delta_sigma_ratio_{ordinal}_MAN",
        ax=ax,
        c=plt.cm.viridis(0.6),
        label="NLLoc",
    )
    ax.axhline(0.683, ls="--", c="k")
    ax.axvline(1, ls="--", c="k")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(ordinal)

    ax.set_xlim(0, 6)
    ax.grid(c="k", lw=1, alpha=0.1)

axs[0].set_ylabel("Proportion of events")
axs[1].set_xlabel("Ratio δ / σ")
axs[2].legend(fontsize=8, loc="lower right")

# Add label
for ax, letter in zip(axs, ["a", "b", "c"]):
    ax.add_artist(
        AnchoredText(
            letter,
            loc="upper right",
            prop={"size": 9, "weight": "bold"},
            frameon=False,
            borderpad=0.1,
        )
    )

fig.tight_layout()

fig.savefig("figureS11.png", dpi=400, bbox_inches="tight")
