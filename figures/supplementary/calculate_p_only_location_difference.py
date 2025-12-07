"""
This script builds Figure S2 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import numpy as np

from obspy.geodetics import gps2dist_azimuth


## Synthetic example
# reference location
lat_ref, lon_ref, z_ref = (-0.000327, -0.001778, 14.95)

# P-only location
lat_p_only, lon_p_only, z_p_only = (-0.000779, -0.001778, 14.75)

dist, ab_az, ba_az = gps2dist_azimuth(lat_ref, lon_ref, lat_p_only, lon_p_only)
z_dist = z_p_only - z_ref

x_dist = dist * np.sin(ab_az * np.pi / 180)
y_dist = dist * np.cos(ab_az * np.pi / 180)

print(f"\nSynthetic location difference: {dist} m horizontal;  {z_dist} m vertical")
print(f"{x_dist} m in X; {y_dist} m in Y\n")


## Askja event example
# reference location
lat_ref, lon_ref, z_ref = (65.016675, -16.655304, 2.5)

# P-only location
lat_p_only, lon_p_only, z_p_only = (65.012944, -16.814355, 2.2)

dist, ab_az, ba_az = gps2dist_azimuth(lat_ref, lon_ref, lat_p_only, lon_p_only)
z_dist = z_p_only - z_ref

x_dist = dist * np.sin(ab_az * np.pi / 180)
y_dist = dist * np.cos(ab_az * np.pi / 180)

print(f"\nAskja event location difference: {dist} m horizontal;  {z_dist} m vertical")
print(f"{x_dist} m in X; {y_dist} m in Y\n")

