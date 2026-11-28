"""int32 scaling constants and overflow checks used when converting strain data to
integer MiniSEED."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import logging

# --- Third party -------------------------------------------------------------
import numpy as np

logger = logging.getLogger(__name__)


INT32_MAX = np.iinfo(np.int32).max

# Scale applied to the data before int32 conversion, per resolution label.
# This table must be the single source for BOTH dfdas_to_mseed (multiplier)
# and the stage_gain written below. Keep only one definition in the module.
INT32_SCALES = {"picostrain": 1e12, "nanostrain": 1e9, "microstrain": 1e6}

# (stage name, output units, units description) for each resolution label
_INT32_STAGE_TEXT = {
    "picostrain": ("Picostrain to Counts Resolution Mapping", "PICOSTRAIN",
                   "Picostrain units (1e-12 strain)"),
    "nanostrain": ("Nanostrain to Counts Resolution Mapping", "NANOSTRAIN",
                  "Nanostrain units (1e-9 strain)"),
    "microstrain": ("Microstrain to Counts Resolution Mapping", "MICROSTRAIN",
                    "Microstrain units (1e-6 strain)"),
}


def check_int32_scaling(data: np.ndarray, scale: float, min_bits: int = 16) -> float:
    """
    Verify that data * scale fits in int32. Raises before any conversion.

    Raises ValueError for NaN/inf, OverflowError if the scaled peak exceeds
    int32. Logs a warning if the scaled peak is so small that the integer
    output has fewer than min_bits of dynamic range.
    Returns the scaled peak.
    """
    arr = np.asarray(data)
    peak = float(np.max(np.abs(arr)))          # NaN/inf propagate through max
    if not np.isfinite(peak):
        n_nan = int(np.isnan(arr).sum())
        n_inf = int(np.isinf(arr).sum())
        raise ValueError(
            f"Cannot convert to int32: data contain {n_nan} NaN and {n_inf} "
            "inf values. Fix or mask them before scaling.")
    scaled_peak = peak * scale
    if scaled_peak > INT32_MAX:
        raise OverflowError(
            f"int32 overflow: peak |value| {peak:.3e} x scale {scale:.0e} = "
            f"{scaled_peak:.3e} exceeds {INT32_MAX}. Largest usable scale is "
            f"{INT32_MAX / peak:.3e}; use a smaller scale/resolution or write "
            "float32 MiniSEED instead.")
    if 0 < scaled_peak < 2 ** min_bits:
        logger.warning(
            "int32 output has little dynamic range: scaled peak is %.1f "
            "counts (< 2^%d). Most of the signal will be lost to rounding; "
            "check the data units and chosen scale.", scaled_peak, min_bits)
    return scaled_peak
