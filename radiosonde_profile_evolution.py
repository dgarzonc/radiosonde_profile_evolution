# ==============================================================================
# vertical_profile_evolution.py
#
# Hovmöller-style diagrams of upper-air radiosonde data over a configurable
# time window. Downloads soundings from the University of Wyoming Upper Air
# Archive, interpolates to a regular 5 hPa pressure grid, and plots the
# temporal evolution of temperature, wind, relative humidity, equivalent
# potential temperature, and their climatological anomalies.
#
# Climatologies are produced by the companion repository and copied manually
# into data/climatology/:
#   https://github.com/dgarzonc/radiosonde_climatology_analysis
#
# Author : David Garzón Casas  —  Physicist, Mg. in Meteorology (Colombia)
# Email  : dgarzonc@unal.edu.co
# ==============================================================================


############
### libs

from datetime import datetime, timedelta
import time
import os
import csv
import traceback
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter
import cmocean.cm as cmo
from metpy import calc as mpcalc
from metpy.units import units
from siphon.simplewebservice.wyoming import WyomingUpperAir
from requests.exceptions import HTTPError, ConnectionError
from socket import gaierror
from urllib3.exceptions import MaxRetryError


############
### paths

# Resolve the repository root as the directory that contains this script,
# so the code works on any machine without manual edits.
path_repo     = os.path.dirname(os.path.abspath(__file__))
folder_data   = f'{path_repo}/data'
folder_output = f'{path_repo}/figures'


# =============================================================================
# Station configuration
#
# Keys are Wyoming platform station codes (used to query Wyoming Upper Air).
# Fields:
#   display_name      : label used in plot titles
#   folder_label      : subfolder name used in data/ and figures/
#                       (must match the climatology folder name copied from
#                        https://github.com/dgarzonc/radiosonde_climatology_analysis)
#   download_hours    : UTC hours for which downloads are attempted
#   valid_hours       : subset of download_hours expected to have data;
#                       used to suppress spurious retry warnings on non-obs times
#   plot_offset_hours : half of the typical inter-sounding interval [h];
#                       used in 'last_days' mode to duplicate the last column
#                       so the most recent sounding has visible width on the plot
# =============================================================================
STATIONS = {
    'SKBO': {
        'display_name':      'Bogotá',
        'folder_label':      'SKBO_Bogota',
        'download_hours':    [12],
        'valid_hours':       [12],
        'plot_offset_hours': 12,   # one sounding per day  → offset = 12 h
    },
    'SKSP': {
        'display_name':      'San Andrés',
        'folder_label':      'SKSP_SanAndres',
        'download_hours':    [0, 12],
        'valid_hours':       [0, 12],
        'plot_offset_hours': 6,    # two soundings per day → offset =  6 h
    },
    'SKLT': {
        'display_name':      'Leticia',
        'folder_label':      'SKLT_Leticia',
        'download_hours':    [12],
        'valid_hours':       [12],
        'plot_offset_hours': 12,
    },
}

# Wyoming codes of stations to process in this run
active_stations = ['SKBO','SKSP','SKLT']

# UTC offset of the local time zone [hours].
# Used to position soundings correctly on the local-time x-axis.
UTC_OFFSET = -5 # For Colombian Local Time: UTC_OFFSET = -5

# Number of most-recent days for which raw files are always re-downloaded
# and re-interpolated on every run, even if they already exist on disk.
# Covers the case where Wyoming updates a sounding with corrected data.
REDOWNLOAD_DAYS = 3


# =============================================================================
# Run mode
#   'last_days'  : process the most recent N_DAYS ending at today's 12 UTC
#   'historical' : process a fixed date range (historical analysis)
# =============================================================================
MODE   = 'last_days' # or 'historical'
N_DAYS = 22           # number of days shown in 'last_days' mode

# Historical analysis parameters — only used when MODE == 'historical'
hist_start_date = datetime(2026, 1, 22, 12, 0)
hist_end_date   = datetime(2026, 2, 18, 12, 0)


############
### logging

def print_(str_to_print):
    """Print a timestamped message to the terminal."""
    line = f"{datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')}: {str_to_print}"
    # _write_log(line)   # log disabled — uncomment to re-enable
    print(line)

# def _write_log(line, mode='a'):
#     log_path = f'{folder_output}/log.txt'
#     try:
#         with open(log_path, mode) as f:
#             if mode == 'a':
#                 line = '\n' + line
#             f.write(line)
#     except Exception as e:
#         print(f'Log write error: {e}')

print_('Script start: vertical_profile_evolution.py')


# =============================================================================
# Resolve date range from run mode
# =============================================================================

if MODE == 'last_days':
    date_end = (datetime.now() + timedelta(hours=5)).replace(
                    hour=12, minute=0, second=0, microsecond=0)
    n_days = N_DAYS

elif MODE == 'historical':
    date_end = hist_end_date.replace(hour=12, minute=0, second=0, microsecond=0)
    n_days   = (hist_end_date.date() - hist_start_date.date()).days + 1

else:
    raise ValueError(f"Unknown MODE '{MODE}'. Use 'last_days' or 'historical'.")

def _format_period_label(date_start, date_end):
    """
    Return an English date-range string for plot titles.

    Same month & year  : "March 1 to 21, 2026"
    Same year only     : "February 12 to March 5, 2026"
    Different years    : "December 25, 2025 to January 15, 2026"
    """
    d0, d1 = date_start, date_end
    if d0.year == d1.year and d0.month == d1.month:
        return f"{d0.strftime('%B')} {d0.day} to {d1.day}, {d1.year}"
    elif d0.year == d1.year:
        return f"{d0.strftime('%B')} {d0.day} to {d1.strftime('%B')} {d1.day}, {d1.year}"
    else:
        return (f"{d0.strftime('%B')} {d0.day}, {d0.year} "
                f"to {d1.strftime('%B')} {d1.day}, {d1.year}")


# Human-readable period label used in plot titles
date_start   = date_end - timedelta(days=n_days - 1)
period_label = _format_period_label(date_start, date_end)

print()
print()
print('=' * 62)
print('  Radiosonde Vertical Profile Evolution')
print('  David Garzon Casas  |  dgarzonc@unal.edu.co')
print('-' * 62)
print(f'  Mode     : {MODE}')
print(f'  Period   : {period_label}')
print(f'  Stations : {", ".join(active_stations)}')
print('=' * 62)
print()
print()


# =============================================================================
# Download
# =============================================================================

############
### load list of known unavailable soundings (avoids re-querying them)

unavailable_path = f'{folder_data}/raw/unavailable_soundings.csv'
os.makedirs(f'{folder_data}/raw', exist_ok=True)
if not os.path.isfile(unavailable_path):
    with open(unavailable_path, 'w', newline='') as f:
        csv.writer(f).writerow(['unavailable_soundings'])
with open(unavailable_path, 'r') as f:
    reader = csv.reader(f)
    unavailable_soundings = [item for row in reader for item in row][1:]

############
### download loop

for sta_code in active_stations:
    sta          = STATIONS[sta_code]
    folder_label = sta['folder_label']
    raw_folder   = f'{folder_data}/raw/{folder_label}'
    updated      = False   # tracks whether any file was (re)downloaded this run

    # Download n_days + 1: the extra day before the period start ensures
    # contourf fills the left edge gap caused by the UTC offset.
    # In historical mode, also download one extra day after date_end
    # to fill the right edge symmetrically.
    dates_to_download = [date_end - timedelta(days=i) for i in range(n_days + 1)]
    if MODE == 'historical':
        dates_to_download += [date_end + timedelta(days=1)]

    for date_i in dates_to_download:

        for hour_i in sta['download_hours']:
            date_i = date_i.replace(hour=hour_i)

            age_days = (datetime.now() - timedelta(hours=UTC_OFFSET) - date_i
                        ).total_seconds() / 86400

            # Always re-download the most recent REDOWNLOAD_DAYS days
            already_on_disk = any(
                date_i.strftime('%Y%m%d_%HZ') in fname
                for fname in os.listdir(raw_folder)
                if fname.endswith('.csv')
            )
            if already_on_disk and (MODE != 'last_days' or age_days > REDOWNLOAD_DAYS):
                continue

            # Skip if previously confirmed as unavailable
            sounding_key = f"{folder_label}_{date_i.strftime('%Y%m%d_%HZ')}"
            if sounding_key in unavailable_soundings:
                continue

            # Attempt download with retries
            max_retries = 15
            retry_delay = 10  # seconds

            for attempt in range(max_retries):
                try:
                    out_fname = f"{date_i.strftime('%Y%m%d_%HZ')}_wyoming.csv"
                    WyomingUpperAir.request_data(date_i, sta_code).to_csv(
                        f'{raw_folder}/{out_fname}', index=False)
                    updated = True
                    break  # success

                except (HTTPError, gaierror, ConnectionError, MaxRetryError):
                    # Silently skip hours that are not expected observation times
                    if hour_i not in sta['valid_hours']:
                        break
                    if attempt < max_retries - 1:
                        print_(f"Server busy — {sta['display_name']} "
                               f"{date_i.strftime('%HZ %d/%m/%Y')}. "
                               f"Retrying in {retry_delay}s "
                               f"({attempt + 1}/{max_retries})")
                        time.sleep(retry_delay)
                    else:
                        print_(f"Failed after {max_retries} attempts: "
                               f"{sta['display_name']} "
                               f"{date_i.strftime('%HZ %d/%m/%Y')}. "
                               f"Try again later.")

                except ValueError:
                    # Sounding not found in Wyoming database
                    if age_days >= 2.0:
                        unavailable_soundings.append(sounding_key)
                    elif hour_i == 12 and date_i.replace(hour=12) == date_end:
                        print_(f"Sounding not found: {sta['display_name']} "
                               f"{date_i.strftime('%HZ %d/%m/%Y')}. "
                               f"Try again later.")
                    break

                except Exception as e:
                    print_(f"Unexpected error: {e}\n{traceback.format_exc()}")
                    break

    if updated:
        print_(f"Data updated: {sta['display_name']} ({sta_code})")

############
### persist updated unavailable-soundings list

unavailable_soundings = (
    ['unavailable_soundings'] + sorted(set(unavailable_soundings), reverse=True)
)
with open(unavailable_path, 'w', newline='') as f:
    writer = csv.writer(f)
    for item in unavailable_soundings:
        writer.writerow([item])

print_('Download complete.')


# =============================================================================
# Interpolation
# Soundings are interpolated onto a regular 5 hPa pressure grid using
# semi-logarithmic interpolation (standard for atmospheric profiles).
# =============================================================================

############
### build pressure grid template

df_template = pd.DataFrame(
    columns=['pressure', 'temperature', 'dewpoint',
             'u_wind', 'v_wind', 'direction', 'speed'])
df_template.loc[len(df_template), 'pressure'] = 1013.25
for p in range(1000, 100 - 5, -5):
    df_template.loc[len(df_template), 'pressure'] = p
df_template['pressure'] = pd.to_numeric(df_template['pressure'], errors='coerce')

############
### interpolation loop

for sta_code in active_stations:
    sta          = STATIONS[sta_code]
    folder_label = sta['folder_label']
    raw_folder   = f'{folder_data}/raw/{folder_label}'
    interp_folder = f'{folder_data}/interpolation/{folder_label}'

    raw_files = [f for f in os.listdir(raw_folder) if f.endswith('.csv')]

    for raw_fname in raw_files:
        date_i        = datetime.strptime(raw_fname[:11], '%Y%m%d_%H')
        interp_fname  = date_i.strftime('%Y%m%d_%HZ.csv')
        interp_path   = f'{interp_folder}/{interp_fname}'

        # Skip if already interpolated, unless in last_days mode and the file
        # is recent enough that it may have been re-downloaded with updated data.
        if os.path.exists(interp_path):
            raw_age_days = (datetime.now() - datetime.fromtimestamp(
                os.path.getmtime(f'{raw_folder}/{raw_fname}'))).total_seconds() / 86400
            if MODE != 'last_days' or raw_age_days > REDOWNLOAD_DAYS:
                continue

        df_raw   = pd.read_csv(f'{raw_folder}/{raw_fname}')
        df_interp = df_template.copy()

        max_press = df_raw['pressure'].max()
        min_press = df_raw['pressure'].min()
        idx_top   = (df_interp['pressure'] - max_press).abs().idxmin()
        idx_bot   = (df_interp['pressure'] - min_press).abs().idxmin()
        press_grid = df_interp['pressure'].loc[idx_top:idx_bot].to_numpy()

        for var in ['temperature', 'dewpoint', 'u_wind', 'v_wind', 'direction', 'speed']:
            press_raw = df_raw['pressure'].to_numpy()
            var_raw   = df_raw[var].to_numpy()
            valid     = ~np.isnan(var_raw)
            press_raw, var_raw = press_raw[valid], var_raw[valid]
            if len(var_raw) < 5:
                continue
            # Semi-logarithmic interpolation (pressure axis in log space)
            var_interp = np.interp(
                np.log(press_grid),
                np.log(press_raw[::-1]),
                var_raw[::-1]
            )
            df_interp.loc[idx_top:idx_bot, var] = var_interp

        df_interp.to_csv(interp_path, index=False)

print_('Interpolation complete.')


# =============================================================================
# Helper functions
# =============================================================================

def _derive_rh(df):
    """Compute relative humidity [%] from dewpoint and temperature.

    Values are clamped to the physical range [0, 100] % after derivation.
    This clamped result is used everywhere RH appears, including anomaly plots.
    """
    sh = mpcalc.specific_humidity_from_dewpoint(
        df['pressure'].to_numpy() * units.hPa,
        df['dewpoint'].to_numpy() * units.degC
    ).to('kg/kg').magnitude
    rh = mpcalc.relative_humidity_from_specific_humidity(
        df['pressure'].to_numpy() * units.hPa,
        df['temperature'].to_numpy() * units.degC,
        sh
    ).to('percent').magnitude
    return np.clip(rh, 0.0, 100.0)


def _derive_ept(df):
    """Compute equivalent potential temperature [°C]."""
    return mpcalc.equivalent_potential_temperature(
        df['pressure'].to_numpy() * units.hPa,
        df['temperature'].to_numpy() * units.degC,
        df['dewpoint'].to_numpy() * units.degC
    ).to('degC').magnitude


def _thin_barbs(df_raw):
    """
    Thin wind barb levels so that no two consecutive levels are within 45 hPa,
    reducing visual clutter on the plot.
    """
    rows_to_delete = []
    row_i = 0
    while row_i in df_raw.index:
        p_ref  = df_raw.loc[row_i, 'pressure']
        nearby = df_raw.loc[
            (p_ref - df_raw['pressure'] < 45) & (p_ref - df_raw['pressure'] > 0)
        ].index.tolist()
        if not nearby:
            row_i += 1
        else:
            rows_to_delete += nearby
            row_i = nearby[-1] + 1
    return df_raw.drop(rows_to_delete).reset_index(drop=True)


def _load_variable(sta_code, date_end, n_days, var, hours, derive_fn=None):
    """
    Build a (pressure × time) data matrix for a given variable.

    Parameters
    ----------
    sta_code  : Wyoming station code (key in STATIONS)
    date_end  : last datetime of the period (at 12 UTC)
    n_days    : number of days to cover
    var       : variable name (column in interpolated CSV, or derived by derive_fn)
    hours     : UTC hours to search for (e.g. [0, 6, 12, 18] or [12])
    derive_fn : optional callable(df) → np.ndarray to compute var from df

    Returns
    -------
    df_data            : DataFrame with 'pressure' + one column per timestamp
    list_dates         : sorted list of datetimes with valid data
    list_dates_nan     : sorted list of datetimes with missing data
    list_days_non_nan  : local dates (at hour=0) that had at least one valid sounding
    """
    folder_label  = STATIONS[sta_code]['folder_label']
    interp_folder = f'{folder_data}/interpolation/{folder_label}'

    # Candidate datetimes: n_days of the period PLUS one extra day before,
    # so that contourf can interpolate backwards and fill the UTC-offset gap
    # at the left edge of the plot.
    # In historical mode, also include one extra day after date_end to fill
    # the right edge symmetrically (last_days uses a ghost column instead).
    days_range = list(range(n_days + 1))             # includes 1 day before start
    extra_after = [date_end + timedelta(days=1)] if MODE == 'historical' else []

    list_dates_all = sorted(
        [(date_end - timedelta(days=d)).replace(hour=h)
         for d in days_range for h in hours]
        + [dt.replace(hour=h) for dt in extra_after for h in hours]
    )

    # Initialize data matrix (pressure × time)
    df_data = pd.DataFrame(
        index=range(len(df_template)),
        columns=['pressure'] + [dt.strftime('%Y%m%d_%HZ') for dt in list_dates_all]
    )
    df_data['pressure'] = df_template['pressure']

    list_days_non_nan = []
    list_dates_nan    = []

    for dt in list_dates_all:
        fname = dt.strftime('%Y%m%d_%HZ.csv')
        if fname in os.listdir(interp_folder):
            df_i   = pd.read_csv(f'{interp_folder}/{fname}')
            values = derive_fn(df_i) if derive_fn else df_i[var].to_numpy()
            df_data[dt.strftime('%Y%m%d_%HZ')] = values
            if int(np.sum(~np.isnan(values.astype(float)))) > 15:
                list_days_non_nan.append((dt + timedelta(hours=UTC_OFFSET)).replace(hour=0))
            else:
                list_dates_nan.append(dt)
        else:
            list_dates_nan.append(dt)

    # Drop columns (soundings) with no data
    df_data = df_data.drop(
        columns=[dt.strftime('%Y%m%d_%HZ') for dt in list_dates_nan],
        errors='ignore'
    )
    list_dates = sorted(set(list_dates_all) - set(list_dates_nan))

    # Drop rows (pressure levels) where all soundings are NaN
    if list_dates:
        df_data = df_data[
            df_data[[dt.strftime('%Y%m%d_%HZ') for dt in list_dates]].notna().any(axis=1)
        ].reset_index(drop=True)

    # In 'last_days' mode: duplicate the most recent sounding column, shifted
    # by plot_offset_hours, so the last day has visible width on the x-axis.
    if MODE == 'last_days' and list_dates:
        offset_h = STATIONS[sta_code]['plot_offset_hours']
        last_dt  = list_dates[-1]
        ghost_dt = last_dt + timedelta(hours=offset_h)
        df_data[ghost_dt.strftime('%Y%m%d_%HZ')] = df_data[last_dt.strftime('%Y%m%d_%HZ')]
        list_dates = list_dates + [ghost_dt]

    return df_data, list_dates, list_dates_nan, list_days_non_nan


def _compute_time_axes(date_end, n_days, list_dates, list_dates_nan, list_days_non_nan):
    """
    Compute normalized day-position vectors for the Hovmöller x-axis.

    The origin (x = 0) is UTC 00:00 of the first displayed day.
    Soundings are positioned in local time (UTC + UTC_OFFSET).
    Missing-day rectangles are centred at local noon.
    """
    origin = date_end.replace(hour=0) - timedelta(days=n_days - 1)

    def _pos_sounding(dt):
        # Convert UTC to local time, then measure distance from origin
        delta = (dt + timedelta(hours=UTC_OFFSET)) - origin
        return (delta.days * 24 + delta.seconds / 3600) / 24

    def _pos_missing(dt):
        delta = dt - origin
        return (delta.days * 24 + delta.seconds / 3600) / 24

    list_days = [_pos_sounding(dt) for dt in list_dates]

    # Identify days with no sounding at all
    seen = set()
    missing_dts = []
    for dt in list_dates_nan:
        local_day = (dt + timedelta(hours=UTC_OFFSET)).replace(hour=0)
        if local_day not in list_days_non_nan and local_day not in seen:
            seen.add(local_day)
            missing_dts.append(local_day + timedelta(hours=12))

    list_no_days = [_pos_missing(dt) for dt in missing_dts]
    return list_days, list_no_days


def _apply_anomaly(df_data, list_dates, sta_code, clim_var, anomaly_type):
    """
    Subtract the climatological daily mean from each sounding column.

    Climatology files are produced by the companion repository
    (https://github.com/dgarzonc/radiosonde_climatology_analysis).
    Each file covers one variable for the full annual cycle:
      filename : {clim_var}.csv   (e.g. temperature.csv)
      columns  : 'pressure', 1, 2, 3, …, 365  (integer day-of-year)
      rows     : pressure levels

    Parameters
    ----------
    clim_var     : variable name that matches the climatology CSV filename,
                   e.g. 'temperature', 'relative_humidity',
                   'equivalent_potential_temperature'
    anomaly_type : 'absolute' → obs − clim  |  'percent' → (obs − clim)/clim × 100

    Returns df_data with anomaly values, or None if the climatology is unavailable.
    """
    folder_label = STATIONS[sta_code]['folder_label']
    clim_folder  = f'{folder_data}/climatology/{folder_label}'

    # Guard: climatology folder must exist
    if not os.path.isdir(clim_folder):
        print_(f"[WARNING] Climatology folder not found for {folder_label}: "
               f"{clim_folder}  —  skipping anomaly plots.")
        return None

    clim_path = f'{clim_folder}/{clim_var}.csv'

    # Guard: climatology file must exist
    if not os.path.isfile(clim_path):
        print_(f"[WARNING] Climatology file not found: {clim_path}  "
               f"—  skipping anomaly plots for {folder_label}.")
        return None

    # Load the full annual-cycle climatology once.
    # Columns after 'pressure' are strings '1'–'365' (pandas reads CSV headers as str).
    df_clim        = pd.read_csv(clim_path)
    clim_pressures = df_clim['pressure'].tolist()

    for dt in list_dates:
        doy     = dt.timetuple().tm_yday   # integer day-of-year (1–365)
        doy_col = str(doy)                 # CSV column header is a string, e.g. '43'

        # Guard: DOY column must exist in the file
        if doy_col not in df_clim.columns:
            print_(f"[WARNING] Day-of-year {doy} not found as a column in "
                   f"{clim_path}  —  skipping anomaly plots for {folder_label}.")
            return None

        col = dt.strftime('%Y%m%d_%HZ')

        for row_i in df_data.index:
            p = df_data.loc[row_i, 'pressure']
            if p in clim_pressures:
                clim_val = df_clim.loc[df_clim['pressure'] == p, doy_col].values[0]
                obs_val  = df_data.loc[row_i, col]
                if anomaly_type == 'absolute':
                    df_data.loc[row_i, col] = obs_val - clim_val
                elif anomaly_type == 'percent':
                    df_data.loc[row_i, col] = (obs_val - clim_val) / clim_val * 100
            else:
                df_data.loc[row_i, col] = np.nan
    return df_data


def _pad_colormap_range(df_data, list_dates, symmetric=True,
                        fixed_vmin=None, fixed_vmax=None):
    """
    Anchor the contourf colormap to the full intended range by inserting two
    phantom pressure levels (p=50 and p=10 hPa) outside the plotted axis,
    carrying the extreme values of the colormap.

    Returns (df_data, vmin, vmax).
    """
    v_absmax = df_data[df_data.columns[1:]].abs().max().max()
    vmax = fixed_vmax if fixed_vmax is not None else v_absmax
    vmin = fixed_vmin if fixed_vmin is not None else (-v_absmax if symmetric else 0)

    n = len(df_data)
    df_data.loc[n,   'pressure'] = 50   # phantom level 1 (empty, just pads y-extent)
    df_data.loc[n+1, 'pressure'] = 10   # phantom level 2 (carries colormap anchors)
    df_data.loc[n+1, list_dates[0].strftime('%Y%m%d_%HZ')] = vmax
    df_data.loc[n+1, list_dates[1].strftime('%Y%m%d_%HZ')] = vmin
    return df_data, vmin, vmax


def _subtitle():
    """Standard credit line for all plots."""
    return (
        f"Issued: {datetime.today().strftime('%d/%m/%Y')}"
        f"  \u2013  Radiosonde data: University of Wyoming"
    )


def _make_title(var_description, sta_code, period_label):
    """Plot title: station and date range on first line, variable on second."""
    return f"{STATIONS[sta_code]['display_name']}  |  {period_label}\n{var_description}"


def _savefig(path):
    """Remove existing file (if any) then save the current figure."""
    if os.path.isfile(path):
        os.remove(path)
    plt.savefig(path, dpi=300, bbox_inches='tight')


def _render_hovmoller(fig, gs, df_data, list_dates, list_days, list_no_days,
                      n_days, date_end, sta_code, title, subtitle,
                      cmap, vmin, vmax, n_levels,
                      label_cbar, norm=None,
                      cbar_top_label=None, cbar_bot_label=None,
                      plot_barbs=True):
    """
    Render the Hovmöller (pressure × time) diagram onto the given figure.

    The plot style — log-pressure y-axis, filled contours + contour lines,
    thinned wind barbs, white rectangles over missing days, and diagonal
    date labels — is preserved exactly from the original script.
    """
    folder_label = STATIONS[sta_code]['folder_label']
    raw_folder   = f'{folder_data}/raw/{folder_label}'

    ax = fig.add_subplot(gs[2:, :])

    # --- Filled contour (Hovmöller)
    contour_kw = dict(levels=n_levels, cmap=cmap, zorder=1)
    if vmin is not None: contour_kw['vmin'] = vmin
    if vmax is not None: contour_kw['vmax'] = vmax
    if norm  is not None: contour_kw['norm'] = norm
    heatmap = ax.contourf(
        list_days, df_data['pressure'],
        df_data[[dt.strftime('%Y%m%d_%HZ') for dt in list_dates]],
        **contour_kw
    )
    ax.contour(heatmap, colors='k', linestyles='solid', alpha=0.2, zorder=2)

    # --- White rectangles over completely missing days
    # Width = 1 day; centred on the missing-day position.
    # Left edge is clamped to 0 (local 00:00 of start date).
    min_press = int(df_data.loc[0, 'pressure'])
    for day_pos in list_no_days:
        left = max(day_pos - 0.5, 0)
        ax.add_patch(patches.Rectangle(
            (left, 100), day_pos + 0.5 - left, min_press - 100,
            linewidth=1, edgecolor='white', facecolor='white', zorder=3
        ))

    # --- Wind barbs (thinned to avoid clutter)
    if plot_barbs:
        for i, dt in enumerate(list_dates):
            raw_path = f"{raw_folder}/{dt.strftime('%Y%m%d_%HZ')}_wyoming.csv"
            if not os.path.exists(raw_path):
                continue
            df_barbs = _thin_barbs(pd.read_csv(raw_path))
            ax.barbs(
                np.ones(len(df_barbs)) * list_days[i],
                df_barbs['pressure'],
                df_barbs['u_wind'], df_barbs['v_wind'],
                length=5.5, pivot='tip', zorder=4
            )

    # --- Log-pressure Y-axis
    ax.invert_yaxis()
    ax.set_yscale('log')
    mandatory_levels = [1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100]
    yticks = [p for p in mandatory_levels if p <= min_press]
    if yticks and abs(yticks[0] - min_press) < 30:
        yticks = yticks[1:]
    yticks = [min_press] + yticks
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(p) for p in yticks])
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.minorticks_off()
    ax.set_ylim([min_press, 100])
    ax.set_ylabel('Pressure [hPa]', fontsize=12)

    # --- X-axis with diagonal date labels
    #
    # All tick marks are placed at local-time midnight (00:00 local) of each
    # displayed day.  The axis origin (x=0) is UTC 00:00 of the first day, so
    # local midnight is offset by  -UTC_OFFSET / 24  axis units (e.g. 5/24 for
    # Colombia UTC-5).  The xlim is set to [local_offset, n_days + local_offset]
    # so the plot area covers exactly local 00:00 of start_date to local 23:59
    # of end_date.
    #
    # For periods > 40 days only 20 evenly-spaced ticks carry a date label
    # (always including start_date and end_date).  Every other daily tick mark
    # is drawn but carries an empty string label.

    local_offset = -UTC_OFFSET / 24.0   # local midnight in axis (UTC) units

    start_date_local = date_end - timedelta(days=n_days - 1)

    # Which day indices (0 = start_date, n_days-1 = end_date) get a text label
    if n_days > 40:
        MAX_TICKS = 20
        raw = np.linspace(0, n_days - 1, MAX_TICKS)
        labeled = sorted(set([0] + [int(round(x)) for x in raw] + [n_days - 1]))
        if len(labeled) > MAX_TICKS:
            interior = [t for t in labeled if t not in (0, n_days - 1)]
            step = max(1, len(interior) // (MAX_TICKS - 2))
            labeled = [0] + interior[::step] + [n_days - 1]
    else:
        labeled = list(range(n_days))

    labeled_set = set(labeled)

    # Build full tick list: every local midnight 0 … n_days (right boundary)
    all_positions = [i + local_offset for i in range(n_days + 1)]
    all_labels    = [
        (start_date_local + timedelta(days=i)).strftime('%d/%m/%Y')
        if i in labeled_set else ''
        for i in range(n_days + 1)   # index n_days = right edge → always ''
    ]

    ax.set_xticks(all_positions)
    ax.set_xticklabels(all_labels, ha='right', va='top', rotation_mode='anchor')
    ax.tick_params(axis='x', rotation=-45)
    ax.set_xticks([], minor=True)
    ax.set_xlim([local_offset, n_days + local_offset])
    for lbl in ax.get_xticklabels():
        lbl.set_horizontalalignment('left')

    # --- Vertical reference line per sounding
    for day_pos in list_days:
        ax.plot([day_pos, day_pos], [100, min_press], c='k', lw=0.4, alpha=0.5)

    # --- Plot borders
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot([0, list_days[-1]], [100, 100], 'k', lw=1.0, zorder=5)

    # --- Colorbar
    cbar = fig.colorbar(heatmap)
    cbar.set_label(label_cbar, fontsize=12)
    if cbar_top_label:
        cbar.ax.text(-0.35, 0.85, cbar_top_label, va='center', ha='center',
                     rotation=90, transform=cbar.ax.transAxes, fontsize=12)
    if cbar_bot_label:
        cbar.ax.text(-0.35, 0.15, cbar_bot_label, va='center', ha='center',
                     rotation=90, transform=cbar.ax.transAxes, fontsize=12)

    # --- Title panel (top strip)
    ax_title = fig.add_subplot(gs[:2, :8])
    ax_title.text(0.5, 0.42, title,
                  ha='center', va='bottom', fontsize=17, fontweight='bold',
                  transform=ax_title.transAxes)
    ax_title.text(0.5, 0.1, subtitle,
                  ha='center', va='bottom', fontsize=12,
                  transform=ax_title.transAxes)
    ax_title.set_xticks([]); ax_title.set_yticks([])
    for spine in ax_title.spines.values():
        spine.set_visible(False)


# =============================================================================
# Plot functions
#
# Each function is a thin wrapper that sets variable-specific parameters
# (colormap, labels, anomaly type, etc.) and delegates to the shared helpers.
# =============================================================================

def plot_temperature(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='temperature', hours=[0, 6, 12, 18])
    if not list_dates: return
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Air Temperature', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.thermal, vmin=None, vmax=None, n_levels=15,
        label_cbar='Air Temperature [\u00b0C]'
    )
    _savefig(f'{folder_out}/1_temperature.jpg')
    plt.close()


def plot_windspeed(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='speed', hours=[0, 6, 12, 18])
    if not list_dates: return
    df_data, vmin, vmax = _pad_colormap_range(
        df_data, list_dates, symmetric=False, fixed_vmin=0)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Wind Speed', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.tempo, vmin=vmin, vmax=vmax, n_levels=15,
        label_cbar='Wind Speed [knots]'
    )
    _savefig(f'{folder_out}/2_speed.jpg')
    plt.close()


def plot_uwind(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='u_wind', hours=[0, 6, 12, 18])
    if not list_dates: return
    df_data, vmin, vmax = _pad_colormap_range(df_data, list_dates, symmetric=True)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Zonal Wind Component (U)', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.curl, vmin=vmin, vmax=vmax, n_levels=15,
        label_cbar     = 'Wind Speed [knots]',
        cbar_top_label = 'Westerly',
        cbar_bot_label = 'Easterly'
    )
    _savefig(f'{folder_out}/3_u_wind.jpg')
    plt.close()


def plot_vwind(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='v_wind', hours=[0, 6, 12, 18])
    if not list_dates: return
    df_data, vmin, vmax = _pad_colormap_range(df_data, list_dates, symmetric=True)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Meridional Wind Component (V)', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.curl, vmin=vmin, vmax=vmax, n_levels=15,
        label_cbar     = 'Wind Speed [knots]',
        cbar_top_label = 'Southerly',
        cbar_bot_label = 'Northerly'
    )
    _savefig(f'{folder_out}/4_v_wind.jpg')
    plt.close()


def plot_relhum(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='rh', hours=[0, 6, 12, 18],
        derive_fn=_derive_rh)
    if not list_dates: return
    df_data, vmin, vmax = _pad_colormap_range(
        df_data, list_dates, symmetric=False, fixed_vmin=0, fixed_vmax=100)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Relative Humidity', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.balance_r, vmin=vmin, vmax=vmax, n_levels=10,
        label_cbar='Relative Humidity [%]'
    )
    _savefig(f'{folder_out}/5_rh.jpg')
    plt.close()


def plot_ept(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='ept', hours=[0, 6, 12, 18],
        derive_fn=_derive_ept)
    if not list_dates: return
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Equivalent Potential Temperature', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.tarn, vmin=None, vmax=None, n_levels=15,
        label_cbar='Equivalent Potential Temperature [\u00b0C]'
    )
    _savefig(f'{folder_out}/6_ept.jpg')
    plt.close()


def plot_anom_temperature(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='temperature', hours=[12])
    if not list_dates: return
    df_data = _apply_anomaly(df_data, list_dates, sta_code,
                              clim_var='temperature', anomaly_type='absolute')
    if df_data is None: return
    df_data, vmin, vmax = _pad_colormap_range(df_data, list_dates, symmetric=True)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Air Temperature Anomaly', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.balance, vmin=vmin, vmax=vmax, n_levels=15,
        norm=colors.CenteredNorm(),
        label_cbar     = 'Absolute Anomaly [\u00b0C]',
        cbar_top_label = 'Warmer than normal',
        cbar_bot_label = 'Cooler than normal'
    )
    _savefig(f'{folder_out}/7_anom_temperature.jpg')
    plt.close()


def plot_anom_rh(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='rh', hours=[12], derive_fn=_derive_rh)
    if not list_dates: return
    df_data = _apply_anomaly(df_data, list_dates, sta_code,
                              clim_var='relative_humidity', anomaly_type='percent')
    if df_data is None: return
    df_data, vmin, vmax = _pad_colormap_range(df_data, list_dates, symmetric=True)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Relative Humidity Anomaly', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.balance_r, vmin=vmin, vmax=vmax, n_levels=15,
        norm=colors.CenteredNorm(),
        label_cbar     = 'Relative Anomaly [%]',
        cbar_top_label = 'More humid than normal',
        cbar_bot_label = 'Drier than normal'
    )
    _savefig(f'{folder_out}/8_anom_rh.jpg')
    plt.close()


def plot_anom_ept(sta_code, date_end, n_days, folder_out, period_label):
    df_data, list_dates, list_dates_nan, list_days_non_nan = _load_variable(
        sta_code, date_end, n_days, var='ept', hours=[12], derive_fn=_derive_ept)
    if not list_dates: return
    df_data = _apply_anomaly(df_data, list_dates, sta_code,
                              clim_var='equivalent_potential_temperature', anomaly_type='percent')
    if df_data is None: return
    df_data, vmin, vmax = _pad_colormap_range(df_data, list_dates, symmetric=True)
    list_days, list_no_days = _compute_time_axes(
        date_end, n_days, list_dates, list_dates_nan, list_days_non_nan)

    fig = plt.figure(figsize=(20, 10))
    gs  = gridspec.GridSpec(nrows=10, ncols=10, figure=fig)
    _render_hovmoller(
        fig, gs, df_data, list_dates, list_days, list_no_days,
        n_days, date_end, sta_code,
        title    = _make_title('Equivalent Potential Temperature Anomaly', sta_code, period_label),
        subtitle = _subtitle(),
        cmap=cmo.balance_r, vmin=vmin, vmax=vmax, n_levels=15,
        norm=colors.CenteredNorm(),
        label_cbar     = 'Relative Anomaly [%]',
        cbar_top_label = 'More unstable than normal',
        cbar_bot_label = 'More stable than normal'
    )
    _savefig(f'{folder_out}/9_anom_ept.jpg')
    plt.close()


# =============================================================================
# Generate plots
# =============================================================================

for sta_code in active_stations:
    folder_label  = STATIONS[sta_code]['folder_label']
    folder_out    = f'{folder_output}/{folder_label}'

    # Create all required subfolders if they do not exist yet
    os.makedirs(f'{folder_data}/raw/{folder_label}',           exist_ok=True)
    os.makedirs(f'{folder_data}/interpolation/{folder_label}', exist_ok=True)
    os.makedirs(f'{folder_data}/climatology/{folder_label}',   exist_ok=True)
    os.makedirs(folder_out,                                     exist_ok=True)

    print_(f"Generating plots for {STATIONS[sta_code]['display_name']} ({sta_code})")

    plot_temperature     (sta_code, date_end, n_days, folder_out, period_label)
    plot_windspeed       (sta_code, date_end, n_days, folder_out, period_label)
    plot_uwind           (sta_code, date_end, n_days, folder_out, period_label)
    plot_vwind           (sta_code, date_end, n_days, folder_out, period_label)
    plot_relhum          (sta_code, date_end, n_days, folder_out, period_label)
    plot_ept             (sta_code, date_end, n_days, folder_out, period_label)
    plot_anom_temperature(sta_code, date_end, n_days, folder_out, period_label)
    plot_anom_rh         (sta_code, date_end, n_days, folder_out, period_label)
    plot_anom_ept        (sta_code, date_end, n_days, folder_out, period_label)

    print_(f"Plots complete for {STATIONS[sta_code]['display_name']}")

print_('Script complete.')
