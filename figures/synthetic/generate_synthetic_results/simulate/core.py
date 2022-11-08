# -*- coding: utf-8 -*-
"""
Small module that provides basic waveform simulations routines.

"""

from typing import Tuple

import numpy as np
from obspy import Trace, Stream, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
import pandas as pd
from quakemigrate.signal.onsets.stalta import sta_lta_centred
from scipy.signal import butter, detrend, lfilter


π = np.pi


class Wavelet:
    """
    Light utility class to encapsulate a Wavelet.

    """

    def __init__(self):
        self.sps
        self.time = None
        self.data = None

    def project(
        self,
        source_polarisation: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Projects the wavelet onto a source polarisation.

        Parameters
        ----------
        source_polarisation: Source polarisation in degrees.

        Returns
        -------
        wavelet_x: Projected x-component of wavelet function.
        wavelet_y: Projected y-component of wavelet function.

        """

        wavelet_x = self.wavelet * np.cos(np.deg2rad(source_polarisation))
        wavelet_y = self.wavelet * np.sin(np.deg2rad(source_polarisation))
        
        return wavelet_x, wavelet_y


class GaussianDerivativeWavelet(Wavelet):

    def __init__(self, frequency, sps, half_timespan):
        """Instantiate GaussianDerivativeWavelet object."""

        delta_T = 1 / frequency
        σ = delta_T / 6
        self.frequency = frequency
        self.sps = sps

        half_timespan += 2*delta_T

        self.time = np.arange(-half_timespan, half_timespan + 1/sps, 1/sps)
        data = -self.time * np.exp(-self.time**2 / (2*σ**2)) / (σ**3 * np.sqrt(2*π))

        # Roll wavelet so first motion is at ~midpoint of array
        self.data = np.roll(data, int(sps*0.5/frequency)) / max(data)


def rotation2d(θ: float) -> np.ndarray:
    """
    Create a 2-D rotation matrix for a given angle, θ (provided in degrees).

    This is defined as:
        | cos(θ) sin(θ)|
        |-sin(θ) cos(θ)|

    Parameters
    ----------
    θ: Angle of rotation of the rotation matrix.

    Returns
    -------
    rotation_matrix: 2 x 2 rotation matrix.

    """

    # Convert from degrees to radians
    θ_r = np.deg2rad(θ)

    return np.array([[np.cos(θ_r), -np.sin(θ_r)],
                     [np.sin(θ_r), np.cos(θ_r)]])


def butter_bandpass(
    lowcut: float,
    highcut: float,
    fs: float,
    order: int = 4
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a butterworth bandpass filter, to be applied using scipy lfilter.

    Parameters
    ----------
    lowcut: Remove frequencies below this value (Hz).
    highcut: Remove frequencies above this value (Hz).
    fs: Sampling rate (Hz).
    order: Filter order (# of corners).
    
    Returns
    -------
    b, a: Numerator & denominator polynomials of the IIR filter.
    
    """
    
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    
    return b, a


def sta_lta_onset(x, filt=[1., 10., 2], onset_params=[0.1, 1.5], fs=50.):
    """
    Calculate the *centred* short-term average / long-term average onset function
    as implemented in QuakeMigrate v1.0.0.as_integer_ratio
                          
                    value assigned here (last sample of LTA window)
                          v
    |---------------------||-------|
               ^               ^
           LTA window      STA window

    Parameters
    ----------
    x : array-like
        Input signal.
    filt : [float, float, int]
        Bandpass filter parameters: [low-cut (Hz), high-cut (Hz), corners]
    onset_params : [float, float]
        [STA window (secs), LTA window (secs)].
    fs : float
        Sampling rate.
    
    Returns
    -------
    onset : array-like
        Onset function.
    filt_x : array-like
        Filtered signal.

    """

    # Detrend
    detrend(x, overwrite_data=True)
    # Filter
    b, a = butter_bandpass(filt[0], filt[1], fs, order=filt[2])
    filt_x = lfilter(b, a, x)
    
    # Calculate onset function
    stw = int(round(onset_params[0]*fs))
    ltw = int(round(onset_params[1]*fs))

    onset = sta_lta_centred(detrend(filt_x), stw, ltw)

    return onset, filt_x


def simulate_waveforms(wavelet, earthquake_coords, lut, magnitude=1, noise=None, angle_of_incidence=0):
    """
    Simulates the waveforms expected for an earthquake within a given LUT.
    Performs simulation in LQT-space (Latitudinal, SV direction, SH direction),
    before rotating onto ZNE based on the ray angles (back-azimuth and
    inclination).

    """

    if noise is None:
        noise = {
            "traveltime": {"P": 0.02, "S": 0.02},
            "amplitude": {"P": 0.1, "S": 0.1}
        }

    inclination = 90 - angle_of_incidence
    earthquake_ijk = lut.index2coord(earthquake_coords, inverse=True)
    
    stream = Stream()
    # Loop over each station and construct the P and S synthetics
    for i, station_data in lut.station_data.iterrows():
        station = station_data["Name"]
        hypo_dist, az, baz = _gps2hypodist_az_baz(
            station_data,
            earthquake_coords,
            lut.unit_conversion_factor
        )
        # amp_factor = 1.
        amp_factor = 10**(magnitude - _attenuate(hypo_dist))

        # Build L component, e.g. the P-phase synthetic
        P = Trace()
        P_ttime = lut.traveltime_to("P", earthquake_ijk, station=station)
        P_ttime += np.random.normal(scale=noise["traveltime"]["P"], size=1)
        roll_by = int(wavelet.sps*P_ttime)
        P_amp_noise = np.random.normal(
            scale=noise["amplitude"]["P"],
            size=len(wavelet.data)
        )
        P.data = np.roll(wavelet.data.copy()*amp_factor*0.5 + P_amp_noise, roll_by)

        # Build Q/T components, e.g. the S-phase synthetic
        S1, S2 = Trace(), Trace()
        S_ttime = lut.traveltime_to("S", earthquake_ijk, station=station)
        S_ttime += np.random.normal(scale=noise["traveltime"]["S"], size=1)
        roll_by = int(wavelet.sps*S_ttime)
        S_amp_noise = np.random.normal(
            scale=noise["amplitude"]["S"],
            size=len(wavelet.data)
        )
        S1.data = np.roll(wavelet.data.copy()*amp_factor + S_amp_noise, roll_by)
        S2.data = np.zeros(len(S1.data)) + S_amp_noise

        lqt_stream = Stream()
        for component, trace in zip("LQT", [P, S1, S2]):
            trace.stats.starttime = UTCDateTime("2021-02-18T12:00:00.0")
            trace.stats.sampling_rate = wavelet.sps
            trace.stats.station = station
            trace.stats.network = "SC"
            trace.stats.channel = f"CH{component}"
            lqt_stream += trace

        zne_stream = lqt_stream.rotate(
            "LQT->ZNE",
            back_azimuth=baz,
            inclination=inclination
        )

        stream += zne_stream

    return stream


def _gps2hypodist_az_baz(
    station_data: pd.DataFrame, earthquake_coords, unit_conversion_factor):
    # Get station location
    stla, stlo, stel = station_data[["Latitude",
                                     "Longitude",
                                     "Elevation"]].values

    # Get event location
    evlo, evla, evdp = earthquake_coords

    # Evaluate epicentral distance between station and event.
    # gps2dist_azimuth returns distances in metres -- magnitudes
    # calculation requires distances in kilometres.
    dist, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
    epi_dist = dist / 1000

    # Evaluate vertical distance between station and event. Convert to
    # kilometres.
    km_cf = 1000 / unit_conversion_factor
    z_dist = (evdp - stel) / km_cf  # NOTE: stel is actually depth.
    
    hypo_dist = np.sqrt(z_dist**2 + epi_dist**2)

    return hypo_dist, az, baz


def _attenuate(dist):
    return 1.11*np.log10(dist/100.) + 0.00189*(dist - 100.) + 3.
