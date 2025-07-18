"""
This script builds Figure 6 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import obspy
import pandas as pd

from quakemigrate.util import DateFormatter


plt.style.use("../../qm_manuscript.mplstyle")
mpl.rcParams["font.family"] = "Helvetica"

# input_data = pathlib.Path.cwd().parents[1] / "data/icequake_example"
input_data = pathlib.Path.cwd() / "generate_results/outputs/runs/paper_run"

# Set QM trigger parameters - MW and MEI are in seconds
marginal_window = 0.1
minimum_event_interval = 0.25
threshold = 3.

# events = pd.read_csv(
#     input_data / "triggered_events/Paper_settings_500Hz_TriggeredEvents.csv"
# )
events = [
    pd.read_csv(infile)
    for infile
    in input_data.glob("trigger/events/paper_run_*_TriggeredEvents*")
]
events = pd.concat(events)
events["CoaTime"] = events["CoaTime"].apply(obspy.UTCDateTime)
events["MinTime"] = events["MinTime"].apply(obspy.UTCDateTime)
events["MaxTime"] = events["MaxTime"].apply(obspy.UTCDateTime)

starttime = obspy.UTCDateTime("2009-020T00:00:00.0")
endtime = starttime + 300.

coalescence = obspy.read(
    str(input_data / "detect/scanmseed/2009_*"),
    starttime=starttime,
    endtime=endtime
)

norm_coa = coalescence.select(station="COA_N")[0]

fig, ax = plt.subplots(1, figsize=(7.08661, 3), constrained_layout=True)

dt = norm_coa.times(type="utcdatetime")
dt = [str(datetime) for datetime in dt]
dt = pd.to_datetime(dt).values
ax.plot(dt, norm_coa.data / 1e5, color="k", lw=0.4, zorder=10)

for _, event in events.iterrows():
    min_dt = event["MinTime"].datetime
    max_dt = event["MaxTime"].datetime
    mw_stt = (event["CoaTime"] - marginal_window).datetime
    mw_end = (event["CoaTime"] + marginal_window).datetime

    if max_dt < starttime.datetime or min_dt > endtime.datetime:
        continue

    ax.axvspan(
        min_dt,
        mw_stt,
        label="Minimum event interval",
        alpha=0.2,
        color="#777777",
        lw=0.001
    )
    ax.axvspan(mw_end, max_dt, alpha=0.2, color="#777777", lw=0.001)
    ax.axvspan(
        mw_stt,
        mw_end,
        label="Marginal window",
        alpha=0.2,
        color="#2ca25f",
        lw=0.001
    )
    ax.axvline(
        event["CoaTime"].datetime,
        label="Triggered event",
        lw=0.25,
        alpha=0.4,
        color="#1F77B4"
    )

ax.axhline(
    threshold,
    label="Detection threshold",
    color="#2c7fb8",
    linestyle="--"
)

ax.set_ylabel("Normalised coalescence")
ax.set_xlabel("Time")
ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S", 2))
ax.set_xlim([starttime.datetime, endtime.datetime])

plt.savefig("figure6.png", dpi=400)
