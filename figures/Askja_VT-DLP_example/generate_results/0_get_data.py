# -*- coding: utf-8 -*-
"""
This script will download the waveform data and an instrument response
inventory from IRIS (in miniSEED and STATIONXML formats, respectively)
for the 24-hour Askja Volcanotectonic (VT) & Deep-Long-Period (DLP) event
example presented in the manuscript:

    QuakeMigrate **

"""

import pathlib

from obspy import UTCDateTime
from obspy.clients.fdsn.mass_downloader import (GlobalDomain, Restrictions,
                                                MassDownloader)

from quakemigrate.io import read_stations

# --- i/o paths ---
station_file = "./inputs/askja_stations.txt"
data_path = pathlib.Path("./inputs/mSEED")

# --- Set time period over which download data ---
starttime = UTCDateTime("2011-298T00:00:00")
endtime = UTCDateTime("2011-301T00:00:00")

# --- Read in station file ---
stations = read_stations(station_file)

# --- Download waveform data and response info
# setup request
restrictions = Restrictions(
    starttime=starttime,
    endtime=endtime,
    chunklength_in_sec=86400,
    network="Z7",
    station=",".join(stations.Name),
    channel_priorities=["HH[ZNE]", "BH[ZNE]"],
    minimum_interstation_distance_in_m=0,
)

# global domain (specifying network and stations instead)
domain = GlobalDomain()

# local archive structure
def get_mseed_storage(network, station, location, channel, starttime,
                      endtime):
    fname = (data_path / f"{starttime.year}" / f"{starttime.julday}" /
        f"{station}_{channel[2]}.m").as_posix()

    return fname

# dataless storage directory
stationxml_storage = "./inputs/DATALESS"

# run download
mdl = MassDownloader(providers=["IRIS"])
mdl.download(domain, restrictions, threads_per_client=4,
             mseed_storage=get_mseed_storage,
             stationxml_storage=stationxml_storage)
