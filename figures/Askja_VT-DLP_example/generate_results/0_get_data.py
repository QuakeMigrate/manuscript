"""
This script will download the waveform data and an instrument response
inventory from IRIS (in miniSEED and STATIONXML formats, respectively)
for the 24-hour Askja volcano (Iceland) Volcanotectonic (VT) & Deep-Long-Period
(DLP) event example presented in the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., nd White, R.S.
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
station_file = "./inputs/askja_stations.txt"
data_path = pathlib.Path("./inputs/mSEED")
stationxml_storage = "./inputs/DATALESS"


# --- Define directory structure for storing waveform data ---
def get_mseed_storage(network, station, location, channel, starttime, endtime):
    fname = (
        data_path
        / f"{starttime.year}"
        / f"{starttime.julday}"
        / f"{station}_{channel[2]}.m"
    ).as_posix()

    return fname


# --- Set network code & client ---
network = "Z7"
datacentres = ["IRIS"]
# global domain (specifying network and stations instead)
domain = GlobalDomain()

# --- Set time period over which download data ---
starttime = UTCDateTime("2011-298T23:50:00")
endtime = UTCDateTime("2011-300T00:10:00")

# --- Read in station file ---
stations = read_stations(station_file)
stations_string = ",".join(stations["Name"])

mdl = MassDownloader(providers=datacentres)

# --- Set up request ---
restrictions = Restrictions(
    starttime=UTCDateTime("2011-298T23:50:00"),
    endtime=UTCDateTime("2011-299T00:00:00"),
    chunklength_in_sec=600,
    network=network,
    station=stations_string,
    channel_priorities=["HH[ZNE]", "BH[ZNE]"],
    minimum_interstation_distance_in_m=0,
)

# --- Download waveform data ---
mdl.download(
    domain,
    restrictions,
    threads_per_client=4,
    mseed_storage=get_mseed_storage,
    stationxml_storage=stationxml_storage,
)

# --- Set up request ---
restrictions = Restrictions(
    starttime=UTCDateTime("2011-299T00:00:00"),
    endtime=UTCDateTime("2011-300T00:00:00"),
    chunklength_in_sec=86400,
    network=network,
    station=stations_string,
    channel_priorities=["HH[ZNE]", "BH[ZNE]"],
    minimum_interstation_distance_in_m=0,
)

# --- Download waveform data ---
mdl.download(
    domain,
    restrictions,
    threads_per_client=4,
    mseed_storage=get_mseed_storage,
    stationxml_storage=stationxml_storage,
)

# --- Set up request ---
restrictions = Restrictions(
    starttime=UTCDateTime("2011-300T00:00:00"),
    endtime=UTCDateTime("2011-300T00:10:00"),
    chunklength_in_sec=600,
    network=network,
    station=stations_string,
    channel_priorities=["HH[ZNE]", "BH[ZNE]"],
    minimum_interstation_distance_in_m=0,
)

# --- Download waveform data ---
mdl.download(
    domain,
    restrictions,
    threads_per_client=4,
    mseed_storage=get_mseed_storage,
    stationxml_storage=stationxml_storage,
)
