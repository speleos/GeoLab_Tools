"""das_processor: processing tools for Distributed Acoustic Sensing (DAS) data."""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("das-processor")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0"

from .config import init_config, load_config

# Public names that live in das_processor.py, loaded lazily so that
# `das-processor init` works without importing obspy, matplotlib, etc.
_LAZY = {
    "DASProcessor": ".das_processor",
    "get_travel_times": ".das_processor",
    "get_coordinates": ".das_processor",
    "get_azimuth_dip": ".das_processor",
    "encode_station_name": ".das_processor",
    "station_xml_creation": ".das_processor",
    "degrade_coordinates": ".das_processor",
    "get_cheby2_coeffs": ".das_processor",
    "zerophase_chebyshev_lowpass_filter": ".das_processor",
    "get_data_from_S3": ".das_processor",
    "baleen_whale_detector": ".das_processor",
}

__all__ = ["init_config", "load_config", "__version__", *_LAZY]


def __getattr__(name):
    if name in _LAZY:
        value = getattr(import_module(_LAZY[name], __name__), name)
        globals()[name] = value  # cache
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_LAZY))
