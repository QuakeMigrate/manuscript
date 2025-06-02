# -*- coding: utf-8 -*-
"""
This script combines the QuakeMigrate and NLLOC outputs for the catalogue of
manually picked earthquakes from the region around Askja volcano (Iceland)
for the purpose of benchmarking the location performance of QuakeMigrate, as
presented in the manuscript:

    QuakeMigrate **

"""

import pathlib

import numpy as np
import pandas as pd
from datetime import timedelta

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth


def read_nlloc_summaryhyp(filename):
    """
    Read key information from NonLinLoc summary hyp file into a
    pandas DataFrame.

    """

    events_df = pd.DataFrame()

    lines = [line.rstrip() for line in open(filename, "r").readlines()]

    lines_start = [i for i, line in enumerate(lines) if line.startswith("NLLOC")]
    lines_end = [i for i, line in enumerate(lines) if line.startswith("END_NLLOC")]

    for start, end in zip(lines_start, lines_end):
        event = read_single_hypocentre(lines[start : end + 1])
        events_df = pd.concat([events_df, event], ignore_index=True)

    events_df["EventID"] = events_df.EventID.astype(int)
    events_df["Err_XYZ"] = np.power(
        (events_df.Err_X * events_df.Err_Y * events_df.Err_Z), 1 / 3
    )

    return events_df


def read_single_hypocentre(lines):
    """
    Read key info from single NLL event in hyp format.

    """

    event = dict()

    lines_dict = dict([line.split(None, 1) for line in lines[:-2]])

    event["FNAME"] = (
        lines_dict["NLLOC"].split()[0].split('"')[1].split("/")[-1] + ".loc.hyp"
    )
    event["LOCATED"] = (
        True if lines_dict["NLLOC"].split()[1].split('"')[1] == "LOCATED" else False
    )
    event["EventID"] = int(lines_dict["PUBLIC_ID"])
    event["OBSFILE"] = lines_dict["SIGNATURE"].split("/")[-1].split()[0]
    try:
        event["DT"] = pd.to_datetime(
            UTCDateTime.strptime(
                " ".join(lines_dict["GEOGRAPHIC"][3:30].split()),
                "%Y %m %d  %H %M %S.%f",
            ).datetime
        )
    except ValueError as e:
        try:
            date = UTCDateTime.strptime(
                " ".join(lines_dict["GEOGRAPHIC"][3:13].split()), "%Y %m %d"
            )
            hours = float(lines_dict["GEOGRAPHIC"][15:17])
            mins = float(lines_dict["GEOGRAPHIC"][18:20])
            secs = float(lines_dict["GEOGRAPHIC"][21:30])
            dt = date + timedelta(hours=hours, minutes=mins, seconds=secs)
            event["DT"] = pd.to_datetime(dt.datetime)
        except ValueError:
            print(e, lines_dict["GEOGRAPHIC"][3:30])
            return
    event["X"] = float(lines_dict["GEOGRAPHIC"].split()[10])
    event["Y"] = float(lines_dict["GEOGRAPHIC"].split()[8])
    event["Z"] = float(lines_dict["GEOGRAPHIC"].split()[12])
    event["RMS"] = float(lines_dict["QUALITY"].split()[7])
    event["NPHA"] = int(lines_dict["QUALITY"].split()[9])
    event["AZ_GAP"] = float(lines_dict["QUALITY"].split()[11])
    event["MIN_DST"] = float(lines_dict["QUALITY"].split()[13])
    event["VPVS"] = float(lines_dict["VPVSRATIO"].split()[1])
    event["NPAIR"] = int(lines_dict["VPVSRATIO"].split()[3])
    event["Err_X"] = np.sqrt(3.53 * float(lines_dict["STATISTICS"].split()[7]))
    event["Err_Y"] = np.sqrt(3.53 * float(lines_dict["STATISTICS"].split()[13]))
    event["Err_Z"] = np.sqrt(3.53 * float(lines_dict["STATISTICS"].split()[17]))

    # read in some info from PHASE lines
    for i, line in enumerate(lines_dict):
        if line == "PHASE":
            phase_line = i
            break

    phase_df = pd.DataFrame([line.split(None) for line in lines[phase_line + 1 :]])
    phase_df.columns = lines[phase_line].split(None)[1:]

    event["NUM_P"] = phase_df.Pha.str.count("P").sum()
    event["NUM_S"] = phase_df.Pha.str.count("S").sum()

    return pd.DataFrame(event, index=[0])


def get_differences(s):
    """Function to apply to each Series to calculate location differences."""

    s["epi_diff"], s["epi_diff_az"], _ = gps2dist_azimuth(
        s.Y_QM, s.X_QM, s.Y_MAN, s.X_MAN
    )
    s["epi_diff"] /= 1000  # convert from m -> km
    s["X_diff"] = s.epi_diff * np.sin(s.epi_diff_az * np.pi / 180)
    s["Y_diff"] = s.epi_diff * np.cos(s.epi_diff_az * np.pi / 180)
    s["Z_diff"] = s.Z_MAN - s.Z_QM

    return s


# main
QM_DIR = pathlib.Path("./outputs/runs/QM-manpick_comparison/locate")
NLL_SUMMARY_HYP = pathlib.Path("./outputs/NLLOC/loc/loc.qm-man_matched.summary.hyp")
OUTFILE = pathlib.Path("./outputs/QM-man_combined.csv")

if __name__ == "__main__":
    # read QM events
    qm_events = pd.DataFrame()
    for evfile in sorted((QM_DIR / "events").glob("*.event")):
        tmp = pd.read_csv(evfile)
        qm_events = pd.concat([qm_events, tmp])
        del tmp
    qm_events.loc[:, "DT"] = pd.to_datetime(qm_events.DT).dt.tz_convert(None)
    qm_events["DT"] = qm_events.DT.astype("datetime64[ns]")
    qm_events = qm_events.set_index("EventID")

    # read NLL events
    nll_events = read_nlloc_summaryhyp(NLL_SUMMARY_HYP)
    nll_events = nll_events.set_index("EventID")

    # combine events
    combined = qm_events.join(
        nll_events[
            [
                "DT",
                "X",
                "Y",
                "Z",
                "RMS",
                "NPHA",
                "NPAIR",
                "NUM_P",
                "NUM_S",
                "MIN_DST",
                "AZ_GAP",
                "VPVS",
                "Err_X",
                "Err_Y",
                "Err_Z",
                "Err_XYZ",
                "FNAME",
                "OBSFILE",
            ]
        ],
        lsuffix="_QM",
        rsuffix="_MAN",
    )

    # calculate location differences
    combined = combined.apply(get_differences, axis=1)

    # save to csv
    combined.to_csv(OUTFILE)

# end
