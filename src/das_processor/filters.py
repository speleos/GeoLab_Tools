"""Zero-phase anti-alias filter design and application."""
from __future__ import annotations

# --- Third party -------------------------------------------------------------
import numpy as np
import scipy.signal as sps
from obspy import Trace


def design_decimation_fir(factor: int, window: str | tuple = ("kaiser", 5.0)) -> np.ndarray:
    """
    The anti-alias FIR that scipy.signal.resample_poly designs by default for
    up=1, down=factor (Kaiser beta=5, 20*factor+1 taps, cutoff at the new Nyquist).
    """
    half_len = 10 * factor
    return sps.firwin(2 * half_len + 1, 1.0 / factor, window=window)


def get_cheby2_coeffs(nyquist: float, freqmax: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculates stable Chebyshev Type II coefficients for the current data.
    Ensures order <= 12 for numerical stability.
    """
    
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    # stop band frequency
    ws = freqmax / nyquist
    # pass band frequency
    wp = ws

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = sps.cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = sps.cheby2(order, rs, wn, btype="low", analog=0, output="ba")
    return b,a


def zerophase_effective_fir(
    b: np.ndarray, a: np.ndarray, tol: float = 1e-6, max_len: int = 200_000
) -> np.ndarray:
    """
    Effective FIR of signal.filtfilt(b, a, x): the autocorrelation of the IIR
    impulse response, truncated where it falls below tol * peak.
    Symmetric, odd length, DC gain ~1.
    """
    imp = np.zeros(max_len)
    imp[0] = 1.0
    h = sps.lfilter(b, a, imp)
    big = np.nonzero(np.abs(h) > tol * np.max(np.abs(h)))[0]
    if big[-1] >= max_len - 1:
        raise ValueError("IIR impulse response did not decay within max_len samples")
    h = h[:big[-1] + 1]
    return sps.fftconvolve(h, h[::-1])


def zerophase_chebyshev_lowpass_filter(trace: Trace, freqmax: float) -> None:
        """
        Custom Chebychev type two zerophase lowpass filter useful for
        decimation filtering.
        This filter is stable up to a reduction in frequency with a factor of
        10. If more reduction is desired, simply decimate in steps.
        Partly based on a filter in ObsPy.
        :param trace: The trace to be filtered.
        :param freqmax: The desired lowpass frequency.
        Will be replaced once ObsPy has a proper decimation filter.
        This code is from LASIF repository (Lion Krischer).
        """
        b, a = get_cheby2_coeffs(trace.stats.sampling_rate * 0.5, freqmax)

        # Apply twice to get rid of the phase distortion.
        trace.data = sps.filtfilt(b, a, trace.data)
