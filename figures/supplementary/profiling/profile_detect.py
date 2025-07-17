"""
This script runs the detect stage for the Askja volcano (Iceland)
Volcanotectonic (VT) & Deep-Long-Period (DLP) event example presented in the
manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

It has been adapted to take the following parameters as arguments on the command line:

    - `n_threads`—the number of threads used for parallel processing
    - `decimation_factor`—the decimation factor used for the LUT
    - `timestep`—the length of time used for each compute stage

"""

import argparse
import pathlib
import sys

from quakemigrate import QuakeScan
from quakemigrate.io import Archive, read_lut, read_stations
from quakemigrate.signal.onsets import STALTAOnset


def main(args: dict | None = None) -> None:
    basepath = pathlib.Path.cwd()
    if not (basepath / "inputs/mSEED").is_dir():
        print(
            "Waveform data for this example have not been downloaded.\n"
            "Please run the Askja VT-DLP scripts in the figures directory first."
        )
        sys.exit(1)

    station_file = basepath / "inputs/askja_stations.txt"
    data_in = basepath / "inputs/mSEED"
    lut_file = basepath / "outputs/lut/askja.LUT"
    run_path = "outputs/runs"
    run_name = "profiling"

    # --- Set time period over which to run detect ---
    starttime = "2011-299T12:00:00.0"
    endtime = "2011-299T13:00:00.0"

    # --- Read in station file ---
    stations = read_stations(station_file)

    # --- Create new Archive and set path structure ---
    archive = Archive(
        archive_path=data_in, stations=stations, archive_format="YEAR/JD/STATION"
    )

    # --- Load the LUT ---
    lut = read_lut(lut_file=lut_file)
    lut.decimate(3 * [args.decimation_factor], inplace=True)

    # --- Create new Onset ---
    onset = STALTAOnset(
        position="classic", sampling_rate=50, signal_transform="env_squared"
    )
    onset.phases = ["P", "S"]
    onset.bandpass_filters = {"P": [2, 16, 2], "S": [2, 14, 2]}
    onset.sta_lta_windows = {"P": [0.2, 1.0], "S": [0.2, 1.0]}

    # --- Create new QuakeScan ---
    scan = QuakeScan(
        archive,
        lut,
        onset=onset,
        run_path=run_path,
        run_name=run_name,
        log=True,
        loglevel="info",
    )

    # --- Set detect parameters ---
    scan.timestep = args.timestep
    scan.threads = args.n_threads

    # --- Run detect ---
    scan.detect(starttime, endtime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--decimation_factor",
        default=2,
        type=int,
    )
    parser.add_argument(
        "--timestep",
        default=300.0,
        type=float,
    )
    parser.add_argument(
        "--n_threads",
        default=16,
        type=int,
    )

    args = parser.parse_args()

    main(args)
