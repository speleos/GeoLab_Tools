"""das_processor: processing tools for Distributed Acoustic Sensing (DAS) data."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("das-processor")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0"

from .config import init_config, load_config

__all__ = ["DASProcessor", "init_config", "load_config", "__version__"]


def __getattr__(name):
    # Lazy import: heavy dependencies load only when DASProcessor is requested
    if name == "DASProcessor":
        from .das_processor import DASProcessor
        return DASProcessor
    raise AttributeError(f"module 'das_processor' has no attribute {name!r}")
