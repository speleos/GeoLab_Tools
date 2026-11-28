"""Station coordinates (cached CSV lookup), cable geometry, station naming and travel times."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import functools
import logging
import os
from typing import Any

# --- Third party -------------------------------------------------------------
import numpy as np
import pandas as pd
from obspy.taup import TauPyModel

# --- Local -------------------------------------------------------------------
from ._log import TRACE

logger = logging.getLogger(__name__)


def get_travel_times(
    event: str, global_params: dict[str, Any], rec_lat: float, rec_lon: float
) -> tuple[list, list]:
    """
    Compute theoretical P and S arrival times for an event at a receiver location.

    Parameters:
    - event: key identifying the source event in the `sources` dictionary.
    - global_params: dict containing at least 'exec_params' and 'sources'.
    - rec_lat, rec_lon: receiver latitude and longitude in degrees.

    Returns:
    - arrivals_p, arrivals_s: lists of TauP arrival objects for P and S phases.
    """
    sources = global_params['sources']
    
    #calculate theoretical arrival time
    source_latitude = sources.get(event, {}).get("latitude")
    source_longitude = sources.get(event, {}).get("longitude")
    source_depth = sources.get(event, {}).get("depth")
    i91 = TauPyModel(model="iasp91")
    arrivals_p = i91.get_travel_times_geo(source_depth_in_km=source_depth, receiver_latitude_in_deg=rec_lat,\
                                              receiver_longitude_in_deg=rec_lon, source_latitude_in_deg=source_latitude,\
                                              source_longitude_in_deg=source_longitude,phase_list=('ttp',), ray_param_tol=0.1)
    arrivals_s = i91.get_travel_times_geo(source_depth_in_km=source_depth, receiver_latitude_in_deg=rec_lat,\
                                              receiver_longitude_in_deg=rec_lon, source_latitude_in_deg=source_latitude,\
                                              source_longitude_in_deg=source_longitude,phase_list=('tts',), ray_param_tol=0.1)
    for arr in arrivals_p:
        arr.time=round(arr.time, 3)

    for arr in arrivals_s:
        arr.time=round(arr.time, 3)
        
    logger.debug("%s %s", arrivals_p, arrivals_s)

    return arrivals_p, arrivals_s


@functools.lru_cache(maxsize=8)
def _load_coords_table(csv_path: str, mtime_ns: int) -> dict[str, dict[Any, dict[str, Any]]]:
    """
    Read the coordinates CSV once and index it by channel.

    Cached on (path, modification time), so editing the file is picked up
    automatically. Returns {"CHANNEL": {channel: row_dict}, "ROI_CHANNEL": {...}};
    only keys whose column exists in the CSV are present. When a channel
    appears more than once the first row wins, as with the old
    df.loc[...].values[0].
    """
    df = pd.read_csv(csv_path)
    tables = {}
    for key in ("CHANNEL", "ROI_CHANNEL"):
        if key in df.columns:
            tables[key] = (df.drop_duplicates(subset=key, keep="first")
                             .set_index(key).to_dict("index"))
    return tables


def clear_coords_cache() -> None:
    """Force the next lookup to re-read the coordinates CSV."""
    _load_coords_table.cache_clear()


def _coords_row(
    channel: int, global_params: dict[str, Any], roi_channel: bool = False
) -> dict[str, Any] | None:
    """Row dict for `channel` from the cached coordinates table, or None."""
    pp = global_params['path_params']
    csv_path = os.path.join(pp['coords_path'], pp['coords_DAS'])
    tables = _load_coords_table(csv_path, os.stat(csv_path).st_mtime_ns)
    key = 'ROI_CHANNEL' if roi_channel else 'CHANNEL'
    if key not in tables:
        raise KeyError(key)
    return tables[key].get(channel)


def get_coordinates(
    channel: int, global_params: dict[str, Any], roi_channel: bool = False, strict: bool = False
) -> tuple[float, float, float]:
    """
    Get coordinates (latitude, longitude, elevation) based on the provided channel.
    The CSV is read once and cached (see _load_coords_table).

    Parameters:
    - channel (int): The channel number for which coordinates are needed.
    - global_params: dictionary with 'exec_params' and 'path_params'
      (coords_path, coords_DAS locate the CSV file).
    - roi_channel (bool): look the channel up in ROI_CHANNEL instead of CHANNEL.
    - strict (bool): raise KeyError for an unknown channel instead of returning
      the hardcoded default position (32.64441, -16.94396, 20).

    Returns:
    Tuple (float, float, float): latitude, longitude, elevation.
    """
    pp = global_params['path_params']

    logger.log(TRACE, "Getting coordinates for channel %s from %s%s",
               channel, pp['coords_path'], pp['coords_DAS'])

    row = _coords_row(channel, global_params, roi_channel)

    if row is None:
        if strict:
            raise KeyError(f"Channel {channel} not found in {pp['coords_DAS']}")
        # Default position for unknown channels (kept for backward compatibility)
        return 32.64441, -16.94396, 20

    return row['LATITUDE'], row['LONGITUDE'], round(row['ALTITUDE'], 0)


def get_azimuth_dip(
    channel: int, global_params: dict[str, Any], roi_channel: bool = False
) -> tuple[float, float]:
    """
    Get azimuth and dip based on the provided channel.
    The CSV is read once and cached (see _load_coords_table).

    Parameters:
    - channel (int): The channel number.
    - global_params: dictionary with 'exec_params' and 'path_params'.
    - roi_channel (bool): look the channel up in ROI_CHANNEL instead of CHANNEL.

    Returns:
    Tuple (float, float): azimuth (0-360) and dip. (0, 0) for an unknown channel.
    """
    pp = global_params['path_params']

    logger.log(TRACE, "Getting azimuth and dip for channel %s from %s%s",
               channel, pp['coords_path'], pp['coords_DAS'])

    row = _coords_row(channel, global_params, roi_channel)

    if row is None:
        azi, dip = 0, 0
    else:
        azi, dip = row['AZIMUTH'], row['DIP']

    if azi < 0:
        azi += 360

    return azi, dip


def encode_station_name(prefix: str, distance: float | str, global_params: dict[str, Any]) -> str:
    """
    Encode a station name using a prefix and a numeric distance.

    Parameters:
    - prefix: string prefix for the station code.
    - distance: numeric distance (metres) used to generate the code suffix.
    - global_params: unused; kept so existing calls keep working.

    Returns:
    - code: formatted station code string.
    """

    distance_m=int(distance)
    code_length = 5 - len(prefix)  # Calculate dynamic length
    factor=10 ** (5-code_length)
    converted_distance = distance_m // factor
    code = f"{prefix}{converted_distance:0{code_length}d}"

    logger.debug(f'Encoded station: name: {code} - distance={distance} code_length={code_length} division_factor={factor}')
    return code


def degrade_coordinates(
    lat: float, lon: float, precision: int = 3, error_metres: float = 200
) -> tuple[float, float]:
    """
    Degrades precision and adds a random spatial displacement (error).
    """
    # 1. Calculate the degree offset for Latitude
    # 1 degree lat ≈ 111,111 meters
    lat_offset = error_metres / 111111.0
    
    # 2. Calculate the degree offset for Longitude
    # 1 degree lon ≈ 111,111 * cos(latitude) meters
    lon_offset = error_metres / (111111.0 * np.cos(np.radians(lat)))
    
    # 3. Apply a random shift within the error range
    # This creates a random displacement in any direction up to the error_meters
    random_lat_shift = np.random.uniform(-lat_offset, lat_offset)
    random_lon_shift = np.random.uniform(-lon_offset, lon_offset)
    
    # 4. Apply shift and then round for precision
    degraded_lat = round(lat + random_lat_shift, precision)
    degraded_lon = round(lon + random_lon_shift, precision)
    
    return degraded_lat, degraded_lon
