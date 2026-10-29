"""das_processor: processing tools for Distributed Acoustic Sensing (DAS) data.

Parts of this package are derived from simpleDAS (ASN-Norway/simpleDAS,
GPL-3.0). See the LICENSE file.

Usage
-----
    das-processor init            # copies defaults.py and sources.py to ./
    from das_processor import *   # as before: process_events, DASProcessor, ...
"""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("das-processor")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0"

from .config import init_config, load_config

# Public names defined in das_processor.py. They are loaded lazily so that
# `das-processor init` works without importing obspy, matplotlib, etc.
_LAZY = {name: ".das_processor" for name in (
    # entry point and validation
    "process_events",
    "create_global_params",
    "validate_source_channel_ranges",
    # main class
    "DASProcessor",
    # DAS file reading/writing (simpleDAS-derived)
    "DASDataFrame",
    "find_DAS_files",
    "load_DAS_files",
    "get_data_indexes",
    "get_filemeta",
    "save_to_DAS_file",
    "str2datetime",
    "convert_index_to_rel_time",
    "create_time_axis",
    "unwrap",
    "wrap",
    # helpers
    "get_travel_times",
    "get_coordinates",
    "get_azimuth_dip",
    "encode_station_name",
    "station_xml_creation",
    "degrade_coordinates",
    "get_cheby2_coeffs",
    "zerophase_chebyshev_lowpass_filter",
    "get_data_from_S3",
    "baleen_whale_detector",
)}

# `from das_processor import *` imports exactly these names.
__all__ = ["init_config", "load_config", "__version__", *_LAZY]


def __getattr__(name):
    if name in _LAZY:
        value = getattr(import_module(_LAZY[name], __name__), name)
        globals()[name] = value  # cache so __getattr__ is not called again
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_LAZY))
