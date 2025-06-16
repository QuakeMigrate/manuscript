"""
This script will download the waveform data and an instrument response
inventory from IRIS (in miniSEED and STATIONXML formats, respectively)
for the Rutford cryoseismicity example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib

from obspy import UTCDateTime
from obspy.clients.fdsn.mass_downloader import (
    GlobalDomain,
    Restrictions,
    MassDownloader,
)

from quakemigrate.io import read_stations


# --- i/o paths ---
station_file = "./inputs/rutford_stations.txt"
data_path = pathlib.Path("./inputs/mSEED")
stationxml_storage = "./inputs/DATALESS"


# --- Define directory structure for storing waveform data ---
def get_mseed_storage(network, station, location, channel, starttime, endtime):
    fname = (
        data_path
        / f"{starttime.year}"
        / f"{starttime.julday:03d}"
        / f"{station}_{channel[2]}.m"
    ).as_posix()

    return fname


# --- Set network code & client ---
network = "YG"
datacentres = ["IRIS"]
# global domain (specifying network and stations instead)
domain = GlobalDomain()

# --- Set time period over which download data ---
time_blocks = [
    ("2009-01-19T23:59:30.0", "2009-01-19T23:59:59.999"),
    # ("2009-01-20T00:00:00.0", "2009-01-20T23:59:59.999"),
    # ("2009-01-21T00:00:00.0", "2009-01-21T23:59:59.999"),
    # ("2009-01-22T00:00:00.0", "2009-01-22T00:00:30.0"),
]

# --- Read in station file ---
stations = read_stations(station_file)
stations_string = ",".join(stations["Name"])

for time_block in time_blocks:
    starttime, endtime = [UTCDateTime(time) for time in time_block]
    # --- Set up request ---
    restrictions = Restrictions(
        starttime=starttime,
        endtime=endtime,
        chunklength_in_sec=86400,
        network=network,
        station=stations_string,
        channel_priorities=["GL[123]"],
        minimum_interstation_distance_in_m=0,
        reject_channels_with_gaps=False,
    )

    # --- Download waveform data ---
    mdl = MassDownloader(providers=datacentres)
    mdl.download(
        domain,
        restrictions,
        threads_per_client=3,
        mseed_storage=get_mseed_storage,
        stationxml_storage=stationxml_storage,
    )
