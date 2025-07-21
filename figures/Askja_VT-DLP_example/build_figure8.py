"""
This script builds Figure 8 of the manuscript:

    Winder, T., Bacon, C.A., Smith, J.D., Hudson, T.S., Drew, J., and White, R.S.
    QuakeMigrate: a Python Package for Automatic Earthquake Detection and Location
    Using Waveform Migration and Stacking. (to be submitted to Seismica).

"""

import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy import read, UTCDateTime
from obspy.geodetics import gps2dist_azimuth


plt.style.use("../../qm_manuscript.mplstyle")
plt.rcParams.update({"font.family": "Helvetica"})


events_df = pd.DataFrame()
for f in sorted(
    pathlib.Path.cwd().glob(
        "generate_results/outputs/runs/paper_run/locate/events/*event"
    )
):
    events_df = pd.concat([events_df, pd.read_csv(f)])
    events_df = events_df.reset_index(drop=True)
events_df.loc[:, "DT"] = pd.to_datetime(events_df.DT).dt.tz_convert(None)


def get_aspect(xmin, xmax, ymin, ymax, zmin, zmax, scale=1000):
    """
    Get aspect ratios for map, lat XS & lon XS for given coords.

    Assumes cross-section depths are in KILOMETRES.

    """

    map_ns_metres, _, _ = gps2dist_azimuth(ymin, xmin, ymax, xmin)
    map_ns_deg = ymax - ymin
    map_ew_metres, _, _ = gps2dist_azimuth(ymin, xmin, ymin, xmax)
    map_ew_deg = xmax - xmin

    map_aspect = (map_ns_metres / map_ns_deg) / (map_ew_metres / map_ew_deg)

    lat_aspect = (map_ns_metres / scale) / map_ns_deg
    lon_aspect = map_ew_deg / (map_ew_metres / scale)

    lon_xs_height_ratio = [(ymax - ymin), (zmax - zmin) / lat_aspect]
    lat_xs_width_ratio = [(xmax - xmin), (zmax - zmin) * lon_aspect]

    map_lon_fig_ratio = sum(lon_xs_height_ratio) * map_aspect / (xmax - xmin)

    map_lat_fig_ratio = (ymax - ymin) / sum(lat_xs_width_ratio) * map_aspect

    return (
        map_aspect,
        lat_aspect,
        lon_aspect,
        lat_xs_width_ratio,
        lon_xs_height_ratio,
        map_lon_fig_ratio,
        map_lat_fig_ratio,
    )


def get_scalebar_deg(dist, xmin, xmax, ymin, ymax):
    """Get length of scalebar for given x-axis scaling."""

    y_mean = np.mean([ymin, ymax])

    map_ew_metres, _, _ = gps2dist_azimuth(y_mean, xmin, y_mean, xmax)
    map_ew_deg = xmax - xmin

    deg = dist * 1000 / (map_ew_metres / map_ew_deg)

    return deg


# Cross-section topography
topo_bedrock_ns = pd.read_csv(
    "./XY_FILES/topo_section_lat_vadalda.xygt",
    sep=r"\s+",
    comment="#",
    names=["X", "Y", "d", "Z"],
)

XY_FILES = pathlib.Path.cwd() / "XY_FILES"

# Askja calderas
askja_calderas = []
for f in sorted(XY_FILES.glob("askjacalderas_*.csv")):
    askja_calderas.append(pd.read_csv(f, comment="#", names=["X", "Y"]))

# Glaciers
glaciers = []
for f in sorted(XY_FILES.glob("big_glaciers/*.csv")):
    glaciers.append(pd.read_csv(f, comment="#", names=["X", "Y"]))

# Central volcanoes
central_volc = []
for f in sorted(XY_FILES.glob("centralvolc/*.csv")):
    central_volc.append(pd.read_csv(f, comment="#", names=["X", "Y"]))

# Fissure swarms
fiss_swarms = []
for f in sorted(XY_FILES.glob("fisswarms/*.csv")):
    fiss_swarms.append(pd.read_csv(f, comment="#", names=["X", "Y"]))

# Station coordinates
stations = pd.read_csv("./generate_results/inputs/askja_stations.txt")

# Iceland coastline
iceland = pd.read_csv("./XY_FILES/iceland_outline.xy", sep=" ", names=["X", "Y"])

## Prep
plt.rcParams.update(
    {
        "font.size": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.labelsize": 8,
        "axes.titlesize": 7,
        "date.autoformatter.hour": "%H:%M",
    }
)

# Make copy of input data
df = events_df.copy().sort_values(by="DT").reset_index(drop=True)

# Set plot extent
xmin = -17.07
xmax = -16.07
ymin = 64.94
ymax = 65.262
zmin = -3
zmax = 25

# Get aspect ratios
(
    map_aspect,
    lat_aspect,
    lon_aspect,
    lat_xs_width_ratio,
    lon_xs_height_ratio,
    map_lon_fig_ratio,
    map_lat_fig_ratio,
) = get_aspect(xmin, xmax, ymin, ymax, zmin, zmax)

# Set height ratio for time axes
time_ax_height = 0.4

# Set height ratio for waveform / spectrogram axes
spec_ax_height = 0.8

# Setup figure
fig_width = 180 / 25.4  # Seismica full-width = 180; half-width = 86 mm
fig_height = (
    map_lat_fig_ratio * fig_width * (1 + time_ax_height + spec_ax_height + 0.15)
)

fig = plt.figure(
    figsize=(fig_width, fig_height), dpi=400, facecolor="w"
)

# GridSpec
gs = GridSpec(
    nrows=3,
    ncols=2,
    width_ratios=lat_xs_width_ratio,
    height_ratios=[1, time_ax_height, spec_ax_height],
)

# gridspec -- outer
spec_gs = GridSpecFromSubplotSpec(
    nrows=2,
    ncols=2,
    height_ratios=[2, 0.8],
    width_ratios=[1, 1],
    subplot_spec=gs[2, :],
    hspace=0.5,
)
# gridspec 1 -- within top subplot
gs1 = GridSpecFromSubplotSpec(
    nrows=2,
    ncols=2,
    height_ratios=[1, 1],
    width_ratios=[1, 1],
    subplot_spec=spec_gs[0, :],
    hspace=0,
)

# Setup axes
ax_map = fig.add_subplot(gs[0, 0], aspect=map_aspect, adjustable="box", anchor="SE")
ax_xs = fig.add_subplot(
    gs[0, 1], aspect=lat_aspect, adjustable="box", anchor="SW", sharey=ax_map
)
ax_time = fig.add_subplot(gs[1, :], anchor="N")

# What to colour by
c = "Z"
vmin = -1
vmax = 22
alpha = 0.8
cmap = sns.color_palette("icefire", as_cmap=True)

cmap_sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))

# Plot map
cmpble = ax_map.scatter(
    df.X,
    df.Y,
    c=df[c],
    cmap=cmap,
    alpha=alpha,
    vmin=vmin,
    vmax=vmax,
    s=((df.ML + 2) ** 5) / 3,
    label="Earthquakes",
)
# Grid
ax_map.grid("k", alpha=0.1)

# Askja nested calderas
for xy_df in askja_calderas:
    ax_map.plot(xy_df.X, xy_df.Y, ls="-", lw=0.5, c="k")

# Stations
ax_map.scatter(
    stations.Longitude,
    stations.Latitude,
    marker="^",
    s=20,
    alpha=0.4,
    color="darkgrey",
    edgecolors="darkgrey",
    label="Stations",
)
ax_map.scatter(
    stations.query("Name == 'SVAD'").Longitude,
    stations.query("Name == 'SVAD'").Latitude,
    marker="^",
    s=80,
    color="limegreen",
    edgecolors="k",
    label="Z7.SVAD",
)

# scalebar
scalebar_len = 10  # Length in kilometres
scalebar = AnchoredSizeBar(
    ax_map.transData,
    size=get_scalebar_deg(scalebar_len, xmin, xmax, ymin, ymax),
    label=f"{scalebar_len} km",
    loc="lower left",
    pad=0.5,
    sep=5,
    frameon=False,
    color="k",
)
ax_map.add_artist(scalebar)


# Colorbar
cax = inset_axes(
    ax_map,
    width="90%",
    height="10%",
    bbox_to_anchor=(0.6, 0.085, 0.4, 0.4),
    bbox_transform=ax_map.transAxes,
    loc="lower left",
)
fig.colorbar(cmap_sm, cax=cax, orientation="horizontal", extend="both")
cax.set_xlabel("Depth b.s.l. / km")
cax.xaxis.set_label_position("top")

# Set axes limits & labels
ax_map.set_xlim(xmin, xmax)
ax_map.set_ylim(ymin, ymax)
ax_map.set_ylabel("Latitude / $\\degree$N")
ax_map.set_xlabel("Longitude / $\\degree$E")
ax_map.tick_params(which="both", bottom=True, top=True, right=True, labelbottom=True)

# Add label
ax_map.add_artist(
    AnchoredText(
        "a",
        loc="upper right",
        prop={"size": 9, "weight": "bold"},
        frameon=False,
        borderpad=0.1,
    )
)

## Plot inset map
inset_xmin = -25
inset_xmax = -13
inset_ymin = 63.25
inset_ymax = 66.7

# Get aspect ratios & make axes
inset_aspect, *_ = get_aspect(xmin, xmax, ymin, ymax, 0, 0)

ax_inset = inset_axes(
    ax_map,
    width="30%",
    height="30%",
    bbox_transform=ax_map.transAxes,
    bbox_to_anchor=(0, 0, 1, 1),
    loc="upper left",
    axes_kwargs={"anchor": "NW", "aspect": inset_aspect, "facecolor": "lightblue"},
    borderpad=0,
)

# Iceland landmass
ax_inset.fill(
    iceland.X,
    iceland.Y,
    facecolor="lightgrey",
    linewidth=0.3,
    linestyle="-",
    edgecolor="k",
)
# Fissure swarms
for xy_df in fiss_swarms:
    ax_inset.fill(xy_df.X, xy_df.Y, facecolor="orange", alpha=0.8)
# Glaciers
for xy_df in glaciers:
    ax_inset.fill(xy_df.X, xy_df.Y, facecolor="ghostwhite")
# Central volcanoes
for xy_df in central_volc:
    ax_inset.plot(xy_df.X, xy_df.Y, ls="-", lw=0.5, c="k")

# Outline of zoom map
ax_inset.plot(
    [xmin, xmin, xmax, xmax, xmin],
    [ymin, ymax, ymax, ymin, ymin],
    c="r",
    lw=1,
    zorder=10,
)

# Tick & axes params
ax_inset.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax_inset.set_xlim(inset_xmin, inset_xmax)
ax_inset.set_ylim(inset_ymin, inset_ymax)

## Plot cross-section
ax_xs.scatter(
    df.Z,
    df.Y,
    c=df[c],
    cmap=cmap,
    alpha=alpha,
    vmin=vmin,
    vmax=vmax,
    s=((df.ML + 2) ** 5) / 3,
)
# Grid
ax_xs.grid("k", alpha=0.1)

# Stations
ax_xs.scatter(
    -stations.query("@xmin < Longitude < @xmax").Elevation,
    stations.query("@xmin < Longitude < @xmax").Latitude,
    marker="<",
    s=20,
    alpha=0.4,
    color="darkgrey",
    edgecolors="darkgrey",
)
ax_xs.scatter(
    -stations.query("Name == 'SVAD'").Elevation,
    stations.query("Name == 'SVAD'").Latitude,
    marker="<",
    s=80,
    color="limegreen",
    edgecolors="k",
)
# Topography
ax_xs.plot(topo_bedrock_ns.Z / -1000, topo_bedrock_ns.Y, c="k", lw=0.5)
ax_xs.fill_betweenx(
    topo_bedrock_ns.Y,
    topo_bedrock_ns.Z / -1000,
    zmax,
    color="grey",
    alpha=0.2,
    zorder=0,
)

# Set axes limits and labels
ax_xs.set_xlim(zmin, zmax)
ax_xs.set_xlabel("Depth below sea level / km")
ax_xs.set_ylim(ymin, ymax)
ax_xs.tick_params(
    which="both", left=True, right=True, top=True, labelleft=False, labelright=False
)

# Add label
ax_xs.add_artist(
    AnchoredText(
        "b",
        loc="upper left",
        prop={"size": 9, "weight": "bold"},
        frameon=False,
        borderpad=0.1,
    )
)

## Legend
handles, labels = ax_map.get_legend_handles_labels()
handles_2, labels_2 = ax_xs.get_legend_handles_labels()
handles.extend(handles_2)
labels.extend(labels_2)

leg = ax_xs.legend(handles, labels, loc="upper right")
leg.set_title("Legend", prop={"weight": "bold"})

## Plot Timeline
ax_time.scatter(
    df.DT,
    df.ML,
    c=df[c],
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    alpha=np.min([alpha * 1.5, 1]),
    s=((df.ML + 2) ** 5) / 3,
    zorder=2,
)
# Grid
ax_time.grid("k", alpha=0.1)

# Set axes labels and limits
ax_time.set_ylabel("M$_L$")
ax_time.set_ylim(-0.7, 1.3)
ax_time.set_xlabel("Time on 2011-10-26 / HH:MM")
ax_time.axhline(0, c="k", lw=1, alpha=0.5)
ax_time.set_xlim(
    UTCDateTime("2011-299T00:00:00").datetime, UTCDateTime("2011-300T00:00:00").datetime
)

# Add # of earthquakes
ax_time.add_artist(
    AnchoredText(f"n = {len(df)}", loc="lower right", borderpad=0.1, frameon=False)
)
# Add label
ax_time.add_artist(
    AnchoredText(
        "c",
        loc="upper left",
        prop={"size": 9, "weight": "bold"},
        frameon=False,
        borderpad=0.1,
    )
)

# Plot cumulative # of events
ax_count = ax_time.twinx()
# Extend to continue to end of plot
ax_count.step(
    np.concatenate(
        (
            [np.datetime64(UTCDateTime("2011-299").datetime)],
            df.DT.values,
            [np.datetime64(UTCDateTime("2011-300").datetime)],
        )
    ),
    np.concatenate(([0], df.index.values + 1, [df.index[-1] + 1])),
    c="k",
    lw=1,
    where="post",
)

# Set axes labels, limits
ax_count.set_ylabel("Cumulative #")
ax_count.set_ylim(bottom=0)

fig.autofmt_xdate(rotation=0, ha="center")

### Waveform & spectrogram plots

qm_out = pathlib.Path("./generate_results/outputs/runs/paper_run/locate")
archive = pathlib.Path("./generate_results/inputs/mSEED")

station = "SVAD"
comp = "Z"
phase = "P"

filt1 = 1.5
filt2 = 30.0  # 60.0
bandpass = True
lowpass = False

pre_time = 3
post_time = 6

FFTs = []
frqs = []

labels = [["d", "f"], ["e", "g"]]

for i, event_id in enumerate(["20111026151311980", "20111026180130180"]):
    pickfile = qm_out / "picks" / f"{event_id}.picks"
    eventfile = qm_out / "events" / f"{event_id}.event"
    picks = pd.read_csv(pickfile)
    event = pd.read_csv(eventfile)
    origin_time = UTCDateTime(event.loc[0, "DT"])
    file = (
        archive
        / f"{origin_time.year}"
        / f"{origin_time.julday:03d}"
        / f"{station}_{comp}.m"
    ).as_posix()
    st = read(file, starttime=origin_time - 20, endtime=origin_time + 60)

    pickline = picks.query("Station == @station and Phase == 'P'").iloc[0]
    if pickline.PickTime == "-1":
        pick = UTCDateTime(pickline.ModelledTime)
    else:
        pick = UTCDateTime(pickline.PickTime)

    s_pickline = picks.query("Station == @station and Phase == 'S'").iloc[0]
    if s_pickline.PickTime == "-1":
        s_pick = UTCDateTime(s_pickline.ModelledTime)
    else:
        s_pick = UTCDateTime(s_pickline.PickTime)

    s_p_time = s_pick - pick

    if len(st) != 1:
        print("Uh oh")
        print(st)
        break
    elif st[0].stats.starttime > pick or st[0].stats.endtime < pick:
        print("Uh oh - no waveform data at time of pick!")
        print(st)
        break

    tr = st[0]

    tr.detrend("demean")
    tr.taper(0.05, "cosine")

    # print "Filtering.."
    if bandpass:
        tr.filter("bandpass", freqmin=filt1, freqmax=filt2, corners=2, zerophase=True)
    if lowpass:
        tr.filter("lowpass", freq=filt2, corners=2, zerophase=True)

    tr_spec = tr.copy().trim(
        starttime=pick - 2 * pre_time, endtime=pick + 2 * post_time
    )
    tr.trim(starttime=pick - pre_time, endtime=pick + post_time)

    sps = tr.stats.sampling_rate
    s_int = 1.0 / sps  # sampling interval
    t = np.arange(0, 1, s_int)  # time vector

    n = tr.stats.npts
    k = np.arange(n)
    secs = k / sps  # for the axis
    T = n / sps
    frq = k / T  # +ve/-ve frequency range
    frq = frq[range(n // 2)]  # just +ve freqency range

    # print "Computing fourier transform..."
    try:
        FFT = (
            np.fft.fft(tr.data) / n
        )  # /n to normalise, since length of trace is max possible frequency
        FFT = abs(FFT[range(n // 2)])  # taking only the +ve half of FFT
        FFT = FFT / np.amax(FFT)  # normalising
    except:
        continue

    FFTs.append(FFT)
    frqs.append(frq)

    ## PLOTTING

    # Add arrow to event on ax_xs
    z = event.Z.values[0]
    y = event.Y.values[0]
    ax_xs.annotate(
        "",
        xy=(z + 0.5, y - 0.0025),
        xytext=(z + 4, y - 0.02),
        arrowprops=dict(arrowstyle="->", color=cmap_sm.to_rgba(z)),
    )

    ax1 = plt.subplot(gs1[0, i])
    ax2 = plt.subplot(gs1[1, i], sharex=ax1)
    ax3 = plt.subplot(spec_gs[1, i])

    # Plot trace
    ax1.plot(
        secs + pre_time,
        tr.data / tr.max(),
        color=cmap_sm.to_rgba(event.Z.values[0]),
        lw=1,
    )
    # Plot picks
    ax1.axvline(2 * pre_time, ls=":", c="r")
    ax1.axvline(2 * pre_time + s_p_time, ls=":", c="dodgerblue")
    ax1.text(
        pre_time * 2 - 0.35, 0.7, "P", transform=ax1.transData, fontsize=8, color="r"
    )
    ax1.text(
        2 * pre_time + s_p_time - 0.35,
        0.7,
        "S",
        transform=ax1.transData,
        fontsize=8,
        color="dodgerblue",
    )
    # Set axes labels / limits
    if i == 0:
        ax1.set_ylabel("Amplitude")
    else:
        ax1.tick_params(labelleft=False)
    ax1.tick_params(right=True, labelbottom=False)
    ax1.set_yticks([-1, 0, 1])
    ax1.set_yticklabels(["", "0", "1"])
    ax1.set_ylim(-1.05, 1.05)
    ax1.grid(c="k", alpha=0.1)
    # Annotate with event stats
    ax1.add_artist(
        AnchoredText(
            f"{labels[i][0]}",
            loc="upper left",
            prop={"size": 9, "weight": "bold"},
            frameon=False,
            borderpad=0.1,
        )
    )
    ax1.add_artist(
        AnchoredText(
            (f"M$_L$ {event.ML.values[0]}\n" f"Z = {event.Z.values[0]} km"),
            loc="lower left",
            prop={"size": 6.5},
            frameon=False,
            borderpad=0.1,
        )
    )

    # Plot spectrogram
    tr_spec.spectrogram(axes=ax2)
    # Plot picks
    ax2.axvline(pre_time * 2, ls=":", c="r")
    ax2.axvline(pre_time * 2 + s_p_time, ls=":", c="dodgerblue")
    # Set axes limits & labels
    ax2.set_xlim(xmin=pre_time, xmax=2 * pre_time + post_time)
    ax2.set_ylim(0, filt2)
    if i == 0:
        ax2.set_ylabel("Freq / Hz")
    else:
        ax2.tick_params(labelleft=False)
    ax2.tick_params(right=True, top=True)
    ax2.set_xlabel("Time / s")
    ax2.set_xticks(np.arange(2 * pre_time - 2, 2 * pre_time + post_time + 1, 2))
    ax2.set_xticklabels(np.arange(-2, 7, 2))
    # Annotate with Trace ID
    ax2.add_artist(
        AnchoredText(
            tr_spec.id,
            loc="upper left",
            frameon=False,
            borderpad=0.1,
            prop={"color": "w", "size": 6.5},
        )
    )

    # Plot FFT
    ax3.plot(frq, FFT / np.max(FFT), color=cmap_sm.to_rgba(event.Z.values[0]), lw=1.0)

    # Axes labels, limits
    ax3.set_xlabel("Frequency / Hz")
    if i == 0:
        ax3.set_ylabel("FFT")
    else:
        ax3.tick_params(labelleft=False)
    ax3.set_xlim(0, filt2)
    ax3.set_ylim(0, 1)
    ax3.grid(c="k", alpha=0.1)
    ax3.add_artist(
        AnchoredText(
            f"{labels[i][1]}",
            loc="upper left",
            prop={"size": 9, "weight": "bold"},
            frameon=False,
            borderpad=0.1,
        )
    )

# ignore warning due to nested subplots
warnings.filterwarnings(action="ignore", category=UserWarning)
fig.tight_layout()

plt.savefig("figure8.png", dpi=400, bbox_inches="tight")
