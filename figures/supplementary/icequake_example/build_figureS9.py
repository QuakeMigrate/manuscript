"""
This script builds Supplementary Figure S9 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


plt.style.use("../../../qm_manuscript.mplstyle")
plt.rcParams.update({"font.family": "Helvetica"})


# Definitions
def get_prec_rec(real, arte, filter_by, filterval, test):
    """Calculate precision and recall for a given filter value."""
    if test == "<":
        tps = len(real[real[filter_by] < filterval])
        fps = len(arte[arte[filter_by] < filterval])
        fns = len(real[real[filter_by] > filterval])
    else:
        tps = len(real[real[filter_by] > filterval])
        fps = len(arte[arte[filter_by] > filterval])
        fns = len(real[real[filter_by] < filterval])
    try:
        prec = tps / (tps + fps)
    except ZeroDivisionError:
        prec = 1.0
    rec = tps / (tps + fns)

    return prec, rec


def plot_prec_rec(real, arte, filter_by, test, ax, labels=False):
    """Plot Precision-Recall curve for labelled data, for a given filter parameter."""

    # Totals
    total_real = len(real[filter_by])
    total_arte = len(arte[filter_by])

    # Ratio
    proportion_real = total_real / (total_real + total_arte)

    # TPR and FPR list initiation
    prec_list = []
    rec_list = []
    filterval_list = []

    # Iterate over all values of filter parameter
    real_min_filterval = real[filter_by].min()
    arte_min_filterval = arte[filter_by].min()
    real_max_filterval = real[filter_by].max()
    arte_max_filterval = arte[filter_by].max()
    filterval_min = min(real_min_filterval, arte_min_filterval)
    filterval_max = max(real_max_filterval, arte_max_filterval)
    for filterval in np.linspace(filterval_min - 0.001, filterval_max + 0.001, 1000):
        prec, rec = get_prec_rec(real, arte, filter_by, filterval, test)
        prec_list.append(prec)
        rec_list.append(rec)
        filterval_list.append(filterval)

    # Calculate AUC
    auc = np.trapezoid(prec_list, rec_list)
    if test == ">":
        auc *= -1
    auc_prec = (auc - proportion_real) / (1 - proportion_real)

    ax.plot(rec_list, prec_list, label=f"AUC = {auc_prec:.4g}", c="royalblue")
    ax.plot(rec_list, np.full(len(rec_list), proportion_real), "--", c="limegreen")

    if labels:
        label_precs = []
        label_recs = []
        for label in labels:
            label_prec, label_rec = get_prec_rec(real, arte, filter_by, label, test)
            label_precs.append(label_prec)
            label_recs.append(label_rec)
            ax.annotate(
                f"{label}",
                xy=(label_rec, label_prec),
                textcoords="offset pixels",
                xytext=(10, 5),
                fontsize=9,
            )
        ax.scatter(label_recs, label_precs, c="k", s=4, label="Sample filter values")

    ax.set_ylim(ymax=1.01)
    ax.set_xlim(-0.01, 1.05)
    ax.set_ylabel("Precision")
    ax.set_xlabel("Recall")
    ax.grid()
    ax.legend(loc="lower left", fontsize=9)

    return auc


# main
labels = pd.read_csv(
    pathlib.Path.cwd()
    / "../../../data/icequake_example/rutford_icequakes_2009-020_0000-0010_labels.csv"
)

fig = plt.figure(figsize=(18 / 2.54, 18 / 2.54), facecolor="w")
axs = fig.subplots(nrows=2)

# plot COV_Err_XYZ panel
plot_prec_rec(
    real=labels.query("Label == 1"),
    arte=labels.query("Label == 0"),
    filter_by="COV_Err_XYZ",
    test="<",
    labels=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5],
    ax=axs[0],
)
# plot COA panel
plot_prec_rec(
    real=labels.query("Label == 1"),
    arte=labels.query("Label == 0"),
    filter_by="COA",
    test=">",
    labels=[10, 9, 8, 7, 6, 5.5, 5, 4.5, 4, 3.5],
    ax=axs[1],
)

# Add label
for ax, letter in zip(axs, ["a", "b"]):
    ax.add_artist(
        AnchoredText(
            letter,
            loc="upper right",
            prop={"size": 9, "weight": "bold"},
            frameon=False,
            borderpad=0.1,
        )
    )

fig.savefig("figureS9.png", dpi=400, bbox_inches="tight")
