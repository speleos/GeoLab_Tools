"""Template-matched baleen whale call detection on ObsPy streams.

Everything here works on a plain ``obspy.Stream``; no ``DASProcessor`` is
needed. The usual entry point is :func:`detect_baleen_whales`::

    from .whale import DetectionParams, detect_baleen_whales

    result = detect_baleen_whales(
        stream, species="fin", detection_mode="auto",
        params=DetectionParams(band_detect=(15, 25)),
        plot=False,
    )
    for d in result.detections:
        print(d["datetime"], d["corr"], d["snr"])

``DASProcessor.baleen_whale_detector`` is a thin wrapper around it.
Detection (:func:`detect_whale_calls`), reporting (:func:`format_report`) and
plotting (:func:`plot_detections`) are separate steps, so each can be used or
tested on its own.
"""
from __future__ import annotations

import datetime
import functools
import logging
import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sps
from matplotlib.figure import Figure
from obspy import Stream, Trace
from obspy.signal.cross_correlation import correlate as obspy_cc

from .config import init_whale_templates

logger = logging.getLogger(__name__)

__all__ = [
    "DetectionParams", "SpectrogramParams", "WhaleResult",
    "validate_stream", "resolve_template", "process_trace",
    "detect_whale_calls", "format_report", "plot_detections",
    "detect_baleen_whales",
]

SPECIES = ("blue", "fin", "humpback", "sei", "minke")
DETECTION_MODES = ("single", "auto", "all")


# ---------------------------------------------------------------------------
# Parameters and result
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DetectionParams:
    """Settings of the envelope-peak + template-correlation detector.

    band_detect:     (low, high) bandpass in Hz applied before detection.
    env_win:         Smoothing window of the amplitude envelope, in seconds.
    thresh:          Amplitude threshold, as a multiple of the median envelope.
    min_dist:        Minimum spacing between candidate peaks, in seconds.
    search_win_sec:  Half-width of the template cross-correlation search, in seconds.
    corr_thresh:     Minimum correlation for a candidate to be kept.
    """

    band_detect: tuple[float, float] = (15, 25)
    env_win: float = 0.2
    thresh: float = 4.5
    min_dist: float = 8.0
    search_win_sec: float = 1.0
    corr_thresh: float = 0.15

    def __post_init__(self) -> None:
        band = tuple(self.band_detect)
        if len(band) != 2 or not 0 < band[0] < band[1]:
            raise ValueError(
                f"band_detect must be (low, high) with 0 < low < high, got {self.band_detect!r}")
        object.__setattr__(self, "band_detect", band)
        for name in ("env_win", "search_win_sec"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)!r}")
        if self.min_dist < 0:
            raise ValueError(f"min_dist must not be negative, got {self.min_dist!r}")


@dataclass(frozen=True)
class SpectrogramParams:
    """Settings of the spectrogram panel in the detection figure.

    source:   "raw" for the unfiltered trace, anything else for the bandpassed one.
    fmin/fmax: Frequency range shown, in Hz.
    nperseg/noverlap/nfft: Spectrogram window length, overlap and FFT size (samples).
    dr:       Dynamic range of the colour scale, in dB.
    cmap:     Matplotlib colormap name.
    """

    source: str = "raw"
    fmin: float = 10
    fmax: float = 40
    nperseg: int = 1042
    noverlap: int = 900
    nfft: int = 2048
    dr: float = 15
    cmap: str = "jet"


@dataclass
class WhaleResult:
    """What :func:`detect_baleen_whales` returns.

    detections: Accepted detections across the processed channel(s), one dict each
                (see :func:`detect_baleen_whales` for the keys).
    processed:  Per-channel working data (filtered trace, envelope, threshold, picks).
    figure:     The detection figure, or None when ``plot`` was off. It is closed
                (no longer managed by pyplot) but can still be saved or re-shown.
    """

    detections: list[dict[str, Any]]
    processed: list[dict[str, Any]]
    figure: Figure | None = None


# ---------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------
def _bandpass(x: np.ndarray, fs: float, lo: float, hi: float, order: int = 4) -> np.ndarray:
    nyq = 0.5 * fs
    lo = max(0.001, float(lo))
    hi = min(float(hi), nyq * 0.99)
    b, a = sps.butter(order, [lo / nyq, hi / nyq], btype="band")
    return sps.filtfilt(b, a, x)


def _rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return np.sqrt(np.mean(x**2))


def _stats_lines(name: str, values: Iterable[float], fmt: str = "{:.3f}") -> list[str]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    labels = ("Mean:    ", "Median:  ", "Std:     ", "min:     ", "max:     ")
    funcs = (np.mean, np.median, np.std, np.min, np.max)
    if values.size == 0:
        return [f"{name} {label}n/a" for label in labels]
    return [f"{name} {label}{fmt.format(func(values))}" for label, func in zip(labels, funcs)]


def _centers_to_edges(
    c: np.ndarray, start: float | None = None, end: float | None = None
) -> np.ndarray:
    c = np.asarray(c, dtype=float)
    if len(c) == 1:
        if start is None:
            start = c[0] - 0.5
        if end is None:
            end = c[0] + 0.5
        return np.array([start, end], dtype=float)

    edges = np.empty(len(c) + 1, dtype=float)
    edges[1:-1] = 0.5 * (c[:-1] + c[1:])
    if start is None:
        start = c[0] - 0.5 * (c[1] - c[0])
    if end is None:
        end = c[-1] + 0.5 * (c[-1] - c[-2])
    edges[0] = start
    edges[-1] = end
    return edges


# ---------------------------------------------------------------------------
# Inputs: stream and template
# ---------------------------------------------------------------------------
def validate_stream(stream: Stream | Trace) -> Stream:
    """Return ``stream`` as a Stream after checking it can be analysed.

    A single Trace is wrapped. Raises TypeError for other types and ValueError
    for an empty stream or one with gaps (masked data).
    """
    if isinstance(stream, Trace):
        stream = Stream([stream])
    if not isinstance(stream, Stream):
        raise TypeError(
            f"stream must be an obspy Stream or Trace, got {type(stream).__name__}")
    if len(stream) == 0:
        raise ValueError("The supplied Stream contains no traces.")
    if any(np.ma.isMaskedArray(tr.data) for tr in stream):
        raise ValueError(
            "The supplied Stream has gaps (masked data). Fill or split them "
            "first, e.g. stream.merge(fill_value=0) or stream.split().")
    return stream


def resolve_template(species: str | None, template: str | None, template_fs: float) -> str:
    """Return the template file to use, from an explicit path or a species name.

    Raises ValueError if neither is given or the species is unknown, and
    FileNotFoundError (after writing the bundled templates to the current
    directory) if the file does not exist.
    """
    if template is None and species is None:
        logger.error("No template or species defined.")
        logger.info("Please check the input parameters and try again.")
        raise ValueError("No template or species defined.")

    if template is None:
        species_lower = species.lower()
        if species_lower in SPECIES:
            template = f"{species_lower.upper()}_{template_fs}Hz.txt"
        else:
            logger.error(f"Unknown species: {species}")
            logger.info("Please check the species name and try again.")
            raise ValueError(f"Unknown species: {species}")

    if not os.path.exists(template):
        logger.error(f"Template file {template} not found")
        logger.info("Please check the file path and try again.")
        init_whale_templates()
        raise FileNotFoundError(
            f"Template file {template} not found. Default templates have been initialized in the current directory.\n"
            "Please check the file path and try again."
        )
    return template


@functools.lru_cache(maxsize=32)
def _load_template(
    template: str, mtime_ns: int, template_fs: float, fs: float, search_win_sec: float
) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    """Read and resample the template once per (file, rates); see _prepare_template."""
    tpl = np.loadtxt(template).astype(np.float64)
    tpl -= tpl.mean()
    if template_fs != fs:
        tpl = sps.resample_poly(
            tpl,
            up=int(round(fs)),
            down=int(round(template_fs))
        ).astype(np.float64)
        tpl -= tpl.mean()

    n_tpl = len(tpl)
    pulse_width = (n_tpl - 1) / fs
    half_pulse_samp = int(round(0.5 * pulse_width * fs))
    n_search_window = int(round(search_win_sec * fs))
    tpl_plus = np.concatenate([np.zeros(n_tpl), tpl, np.zeros(n_tpl)])
    return tpl, tpl_plus, n_tpl, half_pulse_samp, n_search_window


def _prepare_template(
    template: str, template_fs: float, fs: float, search_win_sec: float
) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    """Template arrays for a trace sampled at ``fs``.

    The file is read and resampled once and cached (re-read if the file
    changes), which matters in ``all`` mode where every channel needs it.
    Fresh array copies are returned so callers cannot alter the cache.
    """
    tpl, tpl_plus, n_tpl, half, nwin = _load_template(
        template, os.stat(template).st_mtime_ns, template_fs, fs, search_win_sec)
    return tpl.copy(), tpl_plus.copy(), n_tpl, half, nwin


def _auto_choose_channel(stream: Stream, band_detect: tuple[float, float]) -> int:
    peaks = []
    for tr in stream:
        x_raw = tr.data.astype(np.float64)
        fs = float(tr.stats.sampling_rate)
        x_det = _bandpass(x_raw, fs, *band_detect)
        peaks.append(np.nanmax(np.abs(x_det)))
    return int(np.nanargmax(peaks))


def _select_indices(
    stream: Stream, detection_mode: str, channel_idx: int | None,
    band_detect: tuple[float, float],
) -> list[int]:
    if detection_mode not in DETECTION_MODES:
        raise ValueError("DETECTION_MODE must be 'single', 'auto', or 'all'")
    if detection_mode == "single":
        if channel_idx is None:
            raise ValueError("Set CHANNEL_IDX when DETECTION_MODE = 'single'")
        if channel_idx < 0 or channel_idx >= len(stream):
            raise IndexError(f"CHANNEL_IDX out of range: {channel_idx}")
        return [channel_idx]
    if detection_mode == "auto":
        return [_auto_choose_channel(stream, band_detect)]
    return list(range(len(stream)))


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------
def process_trace(
    tr: Trace, idx: int, template: str, template_fs: float, params: DetectionParams
) -> dict[str, Any]:
    """Detect calls in one trace. Returns its working data and accepted detections."""
    x_raw = tr.data.astype(np.float64)
    fs = float(tr.stats.sampling_rate)
    start_time = tr.stats.starttime.datetime
    station = getattr(tr.stats, "station", "")
    channel_label = station if station else f"channel_{idx}"

    tpl, tpl_plus, n_tpl, half_pulse_samp, n_search_window = \
        _prepare_template(template, template_fs, fs, params.search_win_sec)
    x_det = _bandpass(x_raw, fs, *params.band_detect)
    env = np.abs(x_det)
    win_n = max(1, int(params.env_win * fs))
    amp_smooth = np.convolve(env, np.ones(win_n) / win_n, mode="same")
    thr = params.thresh * np.median(env)
    dist = max(1, int(params.min_dist * fs))
    raw_picks, _ = sps.find_peaks(amp_smooth, height=thr, distance=dist)

    results = []
    for p in raw_picks:
        approx_start = max(0, p - half_pulse_samp)
        i0 = approx_start - n_tpl
        i1 = approx_start + 2 * n_tpl
        if i0 < 0 or i1 > len(x_det):
            continue
        st2 = x_det[i0:i1]
        if len(st2) != len(tpl_plus):
            continue

        cc = obspy_cc(st2, tpl_plus, n_search_window,
                      normalize="naive", method="direct", demean=True)
        i_max = int(np.argmax(cc))
        master_corr = float(np.max(cc))
        master_shift = (i_max - n_search_window) / fs
        if master_corr < params.corr_thresh:
            continue

        refined_start = approx_start + (i_max - n_search_window)
        sig0 = refined_start
        sig1 = refined_start + n_tpl
        noi0 = refined_start - n_tpl
        noi1 = refined_start
        if noi0 < 0 or sig1 > len(x_det):
            continue

        signal_win = x_det[sig0:sig1]
        noise_win = x_det[noi0:noi1]
        amp_max = float(np.max(np.abs(signal_win)))
        signal_rms = float(_rms(signal_win))
        noise_rms = float(_rms(noise_win))
        snr = np.nan if noise_rms == 0 else signal_rms / noise_rms
        tsec = refined_start / fs
        abs_time = start_time + datetime.timedelta(seconds=float(tsec))
        results.append({
            "sample": int(refined_start),
            "time_s": float(tsec),
            "datetime": abs_time,
            "corr": float(master_corr),
            "shift_s": float(master_shift),
            "amp_max": float(amp_max),
            "signal_rms": float(signal_rms),
            "noise_rms": float(noise_rms),
            "snr": float(snr),
        })

    return {
        "idx": idx, "label": channel_label, "fs": fs,
        "start_time": start_time, "x_raw": x_raw, "x_det": x_det,
        "amp_smooth": amp_smooth, "thr": thr,
        "raw_picks": raw_picks, "results": results,
    }


def detect_whale_calls(
    stream: Stream, template: str, template_fs: float = 100,
    detection_mode: str = "single", channel_idx: int | None = None,
    params: DetectionParams | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Run the detector. Pure computation: no printing, plotting or file output.

    ``template`` is a path (see :func:`resolve_template`). Returns
    ``(processed, all_results)``: the per-channel working results, and the flat
    list of accepted detections across channels, each tagged with
    ``channel_idx`` and ``channel_label``.
    """
    params = params or DetectionParams()
    indices = _select_indices(stream, detection_mode, channel_idx, params.band_detect)
    processed = [
        process_trace(stream[idx], idx, template, template_fs, params) for idx in indices
    ]
    all_results = []
    for p in processed:
        for rec in p["results"]:
            rec2 = rec.copy()
            rec2["channel_idx"] = p["idx"]
            rec2["channel_label"] = p["label"]
            all_results.append(rec2)
    return processed, all_results


# ---------------------------------------------------------------------------
# Reporting and plotting
# ---------------------------------------------------------------------------
def _choose_plot_result(processed: list[dict[str, Any]]) -> dict[str, Any]:
    if len(processed) == 1:
        return processed[0]
    kept_counts = [len(p["results"]) for p in processed]
    if np.max(kept_counts) > 0:
        return processed[int(np.argmax(kept_counts))]
    raw_counts = [len(p["raw_picks"]) for p in processed]
    return processed[int(np.argmax(raw_counts))]


def format_report(
    detection_mode: str, processed: list[dict[str, Any]],
    all_results: list[dict[str, Any]], params: DetectionParams,
) -> str:
    """Build the human-readable detection summary (no printing)."""
    band_detect, env_win = params.band_detect, params.env_win
    out: list[str] = []
    if detection_mode in ("single", "auto"):
        p = processed[0]
        out += [
            f"Detection mode       : {detection_mode}",
            f"Chosen channel index : {p['idx']}",
            f"Chosen channel label : {p['label']}",
            f"Sampling rate        : {p['fs']:.1f} Hz",
            f"Frequency BW         : {band_detect} Hz",
            f"ENV_WIN              : {env_win} s",
            f"Raw picks            : {len(p['raw_picks'])}",
            f"Kept detections      : {len(p['results'])}",
        ]
        for i, rec in enumerate(p["results"], start=1):
            out.append(
                f"{i:02d} | t = {rec['time_s']:7.2f} s | {rec['datetime']} | "
                f"Cmax = {rec['corr']:.3f} | Shift = {rec['shift_s']:.3f} s | "
                f"Amp = {rec['amp_max']:.3f} | Signal RMS = {rec['signal_rms']:.3f} | "
                f"Noise RMS = {rec['noise_rms']:.3f} | SNR = {rec['snr']:.2f}"
            )
    else:
        nonzero = [
            p for p in processed
            if len(p["raw_picks"]) > 0 or len(p["results"]) > 0
        ]
        out += [
            f"Detection mode       : {detection_mode}",
            f"Processed channels   : {len(processed)}",
            f"Frequency BW         : {band_detect} Hz",
            f"ENV_WIN              : {env_win} s",
            "",
            f"Channels with picks  : {len(nonzero)}",
            f"Total kept detections: {len(all_results)}",
            "",
        ]
        for p in nonzero:
            out.append(
                f"idx={p['idx']:5d} | label={p['label']} | "
                f"raw={len(p['raw_picks']):3d} | kept={len(p['results']):3d}"
            )

    out += ["*******************************", f"Number of lines found {len(all_results)}"]
    if all_results:
        out += _stats_lines("Cmax", [r["corr"] for r in all_results], "{:.3f}")
        out.append("")
        out += _stats_lines("Shift", [r["shift_s"] for r in all_results], "{:.3f}")
        out.append("")
        out += _stats_lines("SNR", [r["snr"] for r in all_results], "{:.3f}")
    return "\n".join(out)


def plot_detections(
    plot_result: dict[str, Any], spectrogram: SpectrogramParams | None = None
) -> Figure:
    """Draw the detection figure for one processed channel and return it.

    Does not save, show or close the figure.
    """
    sp = spectrogram or SpectrogramParams()
    picks = np.array(
        [rec["sample"] for rec in plot_result["results"]], dtype=int
    )
    raw_picks = plot_result["raw_picks"]
    amp_smooth = plot_result["amp_smooth"]
    thr = plot_result["thr"]
    fs = plot_result["fs"]
    x_raw = plot_result["x_raw"]
    x_det = plot_result["x_det"]
    channel_label = plot_result["label"]
    channel_idx = plot_result["idx"]
    start_time = plot_result["start_time"]

    t = np.arange(len(amp_smooth)) / fs
    x_spec = x_raw if sp.source == "raw" else x_det
    f, tt, Sxx = sps.spectrogram(
        x_spec, fs=fs, window="hann", nperseg=sp.nperseg,
        noverlap=sp.noverlap, nfft=sp.nfft, detrend=False,
        scaling="spectrum", mode="magnitude"
    )
    Sxx_db = 20 * np.log10(Sxx + 1e-12)
    sel = (f >= sp.fmin) & (f <= sp.fmax)
    duration = len(x_spec) / fs
    t_edges = _centers_to_edges(tt, start=0.0, end=duration)
    f_sel = f[sel]
    f_edges = _centers_to_edges(f_sel)
    vmax = np.percentile(Sxx_db[sel, :], 99.9)
    vmin = vmax - sp.dr

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(12, 6), sharex=True,
        gridspec_kw={"height_ratios": [1, 2]}
    )
    fig.subplots_adjust(right=0.80, hspace=0.10)
    ax0.plot(t, amp_smooth, lw=1.0, color="black",
             label="Smoothed amplitude")
    ax0.axhline(thr, color="red", ls="--", lw=1.0,
                label="Amplitude threshold")
    if len(raw_picks):
        ax0.plot(
            raw_picks / fs, amp_smooth[raw_picks],
            "o", ms=4, mfc="none", mec="orange", label="Raw picks"
        )
    ax0.set_ylabel("Amplitude")
    ax0.set_title(
        f"Channel {channel_idx} ({channel_label}) | start = {start_time}"
    )
    ax0.set_xlim(0, duration)
    ax0.legend(
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0, frameon=False
    )
    ax1.pcolormesh(
        t_edges, f_edges, Sxx_db[sel, :],
        shading="flat", cmap=sp.cmap, vmin=vmin, vmax=vmax
    )
    ax1.set_ylim(sp.fmin, sp.fmax)
    ax1.set_xlim(0, duration)
    ax1.set_ylabel("Frequency (Hz)")
    ax1.set_xlabel("Time (s)")
    for p in picks:
        ax1.axvline(p / fs, color="white", ls="--", lw=1.0, alpha=0.9)
    plt.setp(ax0.get_xticklabels(), visible=False)
    return fig


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def detect_baleen_whales(
    stream: Stream | Trace, *,
    species: str | None = None, template: str | None = None, template_fs: float = 100,
    detection_mode: str = "single", channel_idx: int | None = None,
    params: DetectionParams | None = None,
    spectrogram: SpectrogramParams | None = None,
    report: bool = True, plot: bool = True,
    save_path: str | os.PathLike | None = "Whale_detector.png", show: bool = True,
) -> WhaleResult:
    """Detect template-matched baleen whale calls in ``stream``.

    The detector selects one channel or processes all of them, finds candidate
    peaks in a smoothed amplitude envelope, and confirms each by
    cross-correlating it with a whale-call template. The stream is analysed
    as-is and is not modified; traces may have different sampling rates.

    Args:
        stream: ObsPy Stream (or a single Trace) without gaps (masked data).
        species: "blue", "fin", "humpback", "sei" or "minke"; selects the bundled
                 template ``<SPECIES>_<template_fs>Hz.txt`` when ``template`` is None.
        template: Path to a template file (one amplitude per line).
        template_fs: Sampling rate of the template file, in Hz.
        detection_mode: "single" (use ``channel_idx``), "auto" (the channel with
                        the largest filtered peak) or "all".
        channel_idx: Zero-based channel for "single" mode.
        params: Detector settings (default :class:`DetectionParams`).
        spectrogram: Figure settings (default :class:`SpectrogramParams`).
        report: Print the detection summary to stdout.
        plot: Draw the figure for one processed channel.
        save_path: File the figure is saved to; None to skip saving.
        show: Call ``plt.show()``. The figure is closed afterwards either way.

    Returns:
        A :class:`WhaleResult`. Each item of ``detections`` has: ``sample``,
        ``time_s`` (seconds from trace start), ``datetime``, ``corr`` (template
        correlation), ``shift_s``, ``amp_max``, ``signal_rms``, ``noise_rms``,
        ``snr`` (NaN when the noise RMS is zero), ``channel_idx`` and
        ``channel_label`` (the station, or ``channel_<index>`` if empty).

    Raises:
        TypeError / ValueError: invalid stream (wrong type, empty, gappy), no
            template or species, unknown species or detection mode, or "single"
            mode without ``channel_idx``.
        FileNotFoundError: the template file does not exist.
        IndexError: ``channel_idx`` is outside the stream.
    """
    params = params or DetectionParams()
    st = validate_stream(stream)
    template = resolve_template(species, template, template_fs)

    processed, all_results = detect_whale_calls(
        st, template, template_fs, detection_mode, channel_idx, params)

    if report:
        print(format_report(detection_mode, processed, all_results, params))

    fig = None
    if plot:
        fig = plot_detections(_choose_plot_result(processed), spectrogram)
        if save_path is not None:
            fig.savefig(save_path, bbox_inches="tight", pad_inches=0)
        if show:
            plt.show()
        plt.close(fig)
    return WhaleResult(all_results, processed, fig)
