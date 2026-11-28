"""Plotting methods of ``DASProcessor``, provided as a mixin.

``DASProcessor`` inherits from :class:`PlottingMixin`, so ``proc.plot_seismic_section(...)`` and the
other ``plot_*`` methods work as before. They rely on the processor's attributes (``dfdas``,
``exec_params`` ...), so they are not meant to be used on their own.
"""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import datetime
import logging
import os
import subprocess
import warnings

# --- Third party -------------------------------------------------------------
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sps
from matplotlib.collections import LineCollection
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.dates import DateFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from obspy import read, Stream, UTCDateTime
from obspy.geodetics import degrees2kilometers, gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel

# --- Local -------------------------------------------------------------------
from .coords import get_travel_times

logger = logging.getLogger(__name__)


# Colormaps whose midpoint is meaningful (zero for signed data).
_DIVERGING_CMAPS = frozenset({
    'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn',
    'Spectral', 'coolwarm', 'bwr', 'seismic',
    'berlin', 'managua', 'vanimo',          # matplotlib >= 3.10
})

# Colormaps that run from low to high with no meaningful centre.
# Cyclic maps (twilight, hsv) are listed here: their wrap-around point is
# not a "zero", so symmetric limits don't help.
_NON_DIVERGING_CMAPS = frozenset({
    'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo',
    'gray', 'binary', 'Greys', 'Blues', 'Greens', 'Oranges', 'Reds',
    'Purples', 'YlOrRd', 'YlGnBu', 'hot', 'afmhot', 'gist_heat', 'copper',
    'bone', 'cubehelix', 'jet', 'rainbow', 'gist_rainbow', 'nipy_spectral',
    'twilight', 'twilight_shifted', 'hsv',
})


def _is_diverging_cmap(cmap: str | Colormap) -> bool | None:
    """
    True if cmap is a known diverging map, False if a known non-diverging
    map, None if it is in neither set. Handles reversed ('_r') names and
    Colormap objects.
    """
    name = cmap if isinstance(cmap, str) else getattr(cmap, "name", "")
    base = name[:-2] if name.endswith("_r") else name
    if base in _DIVERGING_CMAPS:
        return True
    if base in _NON_DIVERGING_CMAPS:
        return False
    return None


def _auto_clim(
    data: np.ndarray, cmap: str | Colormap, pct: float = 98, symmetric: bool | None = None
) -> tuple[float, float]:
    """
    Colour limits for imshow/pcolormesh that suit the colormap.

    symmetric: True forces (-v, v); False uses the data percentiles;
    None (default) decides from the colormap sets above. A colormap in
    neither set is treated as non-diverging, with a warning.
    """
    if symmetric is None:
        symmetric = _is_diverging_cmap(cmap)
        if symmetric is None:
            warnings.warn(
                f"Colormap {cmap!r} is not classified as diverging or "
                "non-diverging; using data percentiles. Pass symmetric=True "
                "for a diverging map.", stacklevel=2)
            symmetric = False
    if symmetric:
        vmax = np.percentile(np.abs(data), pct)
        return -vmax, vmax
    lo, hi = np.percentile(data, [100 - pct, pct])
    return lo, hi


class PlottingMixin:
    """The ``plot_*`` methods of DASProcessor (see the module docstring)."""

    def plot_channel_map(self, offsets, original_offsets):
        """
        Wraps the offset range into a 2D map, `cols` columns wide, with the
        number of lines calculated automatically. One cell = one offset.
        - Green (1): Requested and present
        - Red (2): Requested but missing
        - Blue (3): Present in original file, but not requested
        - White (0): Nothing at this offset
        """
        if original_offsets is None or len(original_offsets) == 0 or offsets is None or len(offsets) == 0:
            logger.warning("No valid offsets are available for plotting.")
            return

        step = max(1, int(self.roi_step)) # symbol width to improve readability

        # 1. Establish the offset range
        min_offset = min(min(original_offsets), min(offsets)) if len(offsets) > 0 else min(original_offsets)
        max_offset = max(max(original_offsets), max(offsets)) if len(offsets) > 0 else max(original_offsets)
        min_offset = int(min_offset)
        max_offset = int(max_offset)
    
        n_units = (max_offset - min_offset) + step

        # Automatic line width
        if n_units > 10000:
            cols = 1000
        elif n_units >= 5000:
            cols = 500
        else:
            cols = 100
        cols = -(-cols // step) * step             # multiple of step: no symbol straddles two lines
        n_rows = -(-n_units // cols)               # ceil division
        
       
        # 2. Matrix with one cell per offset unit; each offset paints `step` cells
        matrix = np.zeros((n_rows, cols), dtype=int) 

        def paint(idx, code):
            idx = np.asarray(idx, dtype=int)
            rows_ = np.repeat((idx // cols)[:, None], step, axis=1)
            cols_ = (idx % cols)[:, None] + np.arange(step)
            ok = cols_ < cols                      # safety clip at the line end
            matrix[rows_[ok], cols_[ok]] = code

    
        orig_idx = (np.asarray(list(original_offsets), dtype=int) - min_offset) 
        req_idx = (np.asarray(offsets, dtype=int) - min_offset) 

        paint(orig_idx, 3)                         # present, not requested (so far)
        paint(req_idx, 2)                          # requested, missing (so far)
        paint(np.intersect1d(orig_idx, req_idx), 1)  # requested and present
        
        # 3. Wrap into lines of `cols` columns, with blank spacer lines in between
        gap = 1                                    # blank lines between data lines
        stride = 1 + gap
        spaced = np.zeros((n_rows * stride - gap, cols), dtype=int)   # 0 = white
        spaced[::stride] = matrix                  # data lines at 0, 2, 4, ...
             # 4. Plot
        cmap = ListedColormap(['white', 'green', 'red', 'blue'])

        fig, ax = plt.subplots(figsize=(16, max(3, n_rows * 0.15)))
        ax.imshow(spaced, cmap=cmap, vmin=0, vmax=3, aspect='auto',
                  interpolation='nearest')

        
        # x: real offset within the line, y: first offset of each line
        x_step = max(1, cols // 10)
        xticks = np.arange(0, cols, x_step)
        ax.set_xticks(xticks)
        #ax.set_xticklabels([str(int(c)) for c in xticks])
        #ax.set_xlabel(f"Offset within line (0-{cols - 1})")
        
        row_step = max(1, n_rows // 15)
        ytick_rows = np.arange(0, n_rows, row_step)
        ax.set_yticks(ytick_rows * stride)         # account for the spacer lines
        ax.set_yticklabels([str(int(min_offset + r * cols)) for r in ytick_rows])
        ax.set_ylabel("Line start offset")

        ax.legend(handles=[
            Patch(facecolor='green', label='Requested & present'),
            Patch(facecolor='red', label='Requested, missing'),
            Patch(facecolor='blue', label='Present, not requested'),
        ], loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3)

        plt.tight_layout()
        plt.show()

    def plot_cross_correlation_heatmap(self, vmin=None, vmax=None, cmap='RdBu_r'):
        """
        Plots all cross-correlation results as a 2D heatmap/waterfall plot.
        """
        if not hasattr(self, 'cross_correlation_results') or not self.cross_correlation_results:
            logger.warning("No cross-correlation results found. Run cross_correlate() first.")
            return

        # 1. Extract channels and sort them to ensure proper spatial ordering
        channels = sorted(list(self.cross_correlation_results.keys()))

        # 2. Grab lags from the first channel to set up our X-axis extent
        sample_channel = channels[0]
        lags = self.cross_correlation_results[sample_channel]['lags']

        # 3. Build a 2D matrix (Shape: num_channels x num_lags)
        corr_matrix = np.array([self.cross_correlation_results[ch]['correlation'] for ch in channels])

        # 4. Set up spatial axis (convert channel indices to actual meters if dx exists)
        dx = self.dfdas.meta.get('dx', 1) if self.dfdas is not None else 1
        distances = [int(ch * dx) for ch in channels]

        # 5. Plotting
        plt.figure(figsize=(12, 7))

        # extent defines [xmin, xmax, ymin, ymax] for the pixel mapping
        extent = [lags.min(), lags.max(), distances[-1], distances[0]]

        # Dynamic color scaling based on data if not provided
        if vmin is None or vmax is None:
            auto_min, auto_max = _auto_clim(corr_matrix, cmap)
            vmin = auto_min if vmin is None else vmin
            vmax = auto_max if vmax is None else vmax       

        im = plt.imshow(corr_matrix, aspect='auto', extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)

        plt.colorbar(im, label='Cross-Correlation Amplitude')
        plt.title('DAS Cross-Correlation Matrix (Waterfall Plot)', fontsize=14, pad=15)
        plt.xlabel('Lag Time (seconds)', fontsize=12)
        plt.ylabel('Distance (metres)', fontsize=12)
        plt.grid(alpha=0.3, linestyle='--')

        plt.tight_layout()
        plt.show()

    def plot_cross_correlation_lines(self, channels_to_plot=None):
        """
        Plots cross-correlation functions as line graphs for specific channels.
        """
        if not hasattr(self, 'cross_correlation_results') or not self.cross_correlation_results:
            logger.warning("No cross-correlation results found.")
            return

        # If no channels are specified, just grab 3 evenly spaced ones
        all_channels = sorted(list(self.cross_correlation_results.keys()))
        if channels_to_plot is None:
            idx = np.linspace(0, len(all_channels) - 1, 3, dtype=int)
            channels_to_plot = [all_channels[i] for i in idx]

        plt.figure(figsize=(10, 5))

        for ch in channels_to_plot:
            if ch in self.cross_correlation_results:
                lags = self.cross_correlation_results[ch]['lags']
                corr = self.cross_correlation_results[ch]['correlation']

                dx = self.dfdas.meta.get('dx', 1) if self.dfdas is not None else 1
                distance = int(ch * dx)

                plt.plot(lags, corr, label=f"Ch {ch} ({distance} m)")

        plt.title('Cross-Correlation Profiles', fontsize=14)
        plt.xlabel('Lag Time (seconds)', fontsize=12)
        plt.ylabel('Correlation Amplitude', fontsize=12)
        plt.axvline(0, color='black', linestyle=':', alpha=0.5, label='Zero Lag')
        plt.legend(loc='upper right')
        plt.grid(alpha=0.4)

        plt.tight_layout()
        plt.show()

    def plot_spectrogram(self, freq_type="no_filter", suffix=""):
        pp = self.path_params
        sfp = self.exec_params.get('spectrogram_freq_params', {})

        if self.mseed_stream is None:
            logger.warning("No ObsPy Stream found in memory. Attempting to create from self.dfdas...")
            self.dfdas_to_mseed(self.event, write=False, format=None, dType=float, decimate=False)

        trace = self.mseed_stream[0].copy()  # Work with a copy to avoid modifying the original

        if freq_type in sfp:
            fmin = sfp[freq_type]["fmin"]
            fmax = sfp[freq_type]["fmax"]
            if fmax == "nyquist":
                fmax=trace.stats.sampling_rate * 0.5
        else:
            fmin=0.01
            fmax=trace.stats.sampling_rate * 0.5

        if fmax>trace.stats.sampling_rate * 0.5:
            fmax=trace.stats.sampling_rate * 0.5

        fmax=40

        trace.filter('bandpass',freqmin=fmin,freqmax=fmax, corners=4, zerophase=True)    

        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_axes([0.1, 1.1, 0.7, 0.2])  # [left, bottom, width, height]
        if freq_type != "whales":
            cmap='coolwarm'
            ax2 = fig.add_axes([0.1, 0.7, 0.7, 0.3]) #low freq
            ax3 = fig.add_axes([0.1, 0.1, 0.7, 0.6]) #high freq
        else:
            cmap='viridis'
            ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.8]) #low freq

        ax4 = fig.add_axes([0.83, 0.1, 0.03, 0.9])

        # Make time vector
        t = np.arange(trace.stats.npts) / trace.stats.sampling_rate
        t_hours = t / 3600

        # Plot waveform (top subfigure)
        ax1.plot(t, trace.data, 'k')
        ax1.set_xlim(t[0], t[-1])

        # Plot spectrogram - Low frequencies (bottom subfigure)
        spec_high_freq = ax2.specgram(trace.data, Fs=trace.stats.sampling_rate, cmap=cmap, NFFT=256, noverlap=128)
        plt.colorbar(spec_high_freq[3], cax=ax4)

        # Plot spectrogram - High frequencies (middle subfigure)
        if freq_type != "whales":
            #spec_low_freq = ax3.specgram(trace.data, Fs=trace.stats.sampling_rate, cmap=cmap, NFFT=1024, noverlap=512)
            ax3.imshow(spec_high_freq[3].get_array(), aspect='auto', extent=spec_high_freq[3].get_extent(), origin='lower', cmap=cmap)
            plt.colorbar(spec_high_freq[3], cax=ax4)

        # Set y-axis scale to logarithmic

        yticks = [0.01, 0.05, 0.1, 0.2, 0.5, 1]
        #ax3.set_ylim(bottom=0.01)  # Set lower limit for logarithmic scale
        if freq_type != "whales":
            ax3.set_yscale('log')
            ax3.set_ylim(bottom=fmin)  # Set lower limit for logarithmic scale
            ax3.set_ylim(top=1)  # Set upper limit for logarithmic scale
            ax3.set_yticks(yticks)
            ax3.spines['top'].set_visible(False)
            ax3.set_yticklabels(yticks)
            ax3.set_xlabel(f'Time (s)')
        else:
            ax1.set_xlabel('Time (s)')

        ax2.set_xticks([])
        ax2.set_ylim(bottom=1)  # Set lower limit for linear scale
        ax2.set_ylim(top=fmax)  # Set lower limit for linear scale
        ax2.spines['bottom'].set_visible(False)
        # Set the x-axis of ax3 to hours
        ax1.set_xticks([])
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel(f'Frequency (Hz)')
        ax4.set_ylabel('Amplitude (dB)')

        if suffix:
            suffix = f"_{suffix}"
        filename=os.path.join(pp.get("figures_path"), f'{self.event}_fmin{fmin}_single_channel_spectrogram_dist-{trace.stats.distance}_{int((t[-1]+1)/3600)}h{suffix}')
        logger.info("%s %s", "saving figure to ", filename)
        distance_km = trace.stats.distance / 1000
        fig.suptitle(f'Spectrogram for channel at {distance_km:.3f} km', y=1.0)
        fig.savefig(filename+'.pdf',dpi=600, format='pdf', bbox_inches='tight', pad_inches=0.1)
        fig.savefig(filename+'.png',dpi=300, format='png')
        plt.show()
        plt.close()  

        return

    def plot_hilbert_envelope(self, data, filename=None):
        pp= self.path_params
        fig, ax = plt.subplots(figsize=(4, 6))

        ax.imshow(data.T, aspect='auto',
                  extent=[self.dfdas.meta["dx"]*self.roi_step*self.min_channel[0]/1000,
                          self.dfdas.meta["dx"]*self.roi_step*self.max_channel[0]/1000,
                          0,
                          self.dfdas.meta["dt"]*self.dfdas.meta["dimensionSizes"][0]],
                          origin='lower',
                          norm="log",
                          cmap='Greys',
                          vmin=0.7,
                          vmax=2.0)

        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Time (s)')
        ax.set_title(f"{filename}\nStart time: {str(self.dfdas.meta['time'])}")

        plt.tight_layout()

        # 6. Save logic
        save_path = os.path.join(pp.get("figures_path"), "Envelopes")
        if not os.path.exists(save_path): os.makedirs(save_path)
        
        filename = os.path.join(save_path, f"env_{self.event.replace(' ', '_')}")
        fig.savefig(f"{filename}.pdf", dpi=300, format='pdf', bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Envelope plot saved as {filename}.pdf")

    def plot_seismic_2d(self, channel_step=1, start_offset_km=30.8,
                        end_offset_km=50.0, dx_meters=12.5, max_time_s=45.0,
                        station=None, event=None, suffix="", **kwargs):
        """
        2D version of seismic plot using imshow.
        """

        if self.dfdas is None:
            logger.warning("No data in memory to plot.")
            return

        event = event or self.event
        logger.info(f"»»» Generating 2D Record Section for {event}...")

        data = self.dfdas.to_numpy()
        dt = self.dfdas.meta['dt']
        fs = 1.0 / dt
        window_start_utc = UTCDateTime(self.dfdas.meta['time'])
        time_axis = (self.dfdas.index - self.dfdas.index[0]).total_seconds()
        
        dx_native = self.dfdas.meta.get('dx', 1.0213)
        absolute_offsets_km = (self.dfdas.columns.to_numpy() * dx_native) / 1000.0

        spatial_mask = (absolute_offsets_km >= start_offset_km) & (absolute_offsets_km <= end_offset_km)
        sub_data = data[:, spatial_mask]
        offsets_km = absolute_offsets_km[spatial_mask]

        if sub_data.shape[1] == 0:
            logger.error(f"ERROR: no channels in range {start_offset_km}-{end_offset_km} km.")
            logger.info(f"Available range: {absolute_offsets_km.min():.2f} - {absolute_offsets_km.max():.2f} km")
            return

        time_mask = time_axis <= max_time_s
        sub_data = sub_data[time_mask, :]
        plot_time = time_axis[time_mask]
        logger.info(f"-> Data: {sub_data.shape[1]} channels, {sub_data.shape[0]} time samples")

        lowcut, highcut = 0.5, 20.0
        nyq = 0.5 * fs
        if highcut >= nyq:
            highcut = nyq * 0.95
        b, a = sps.butter(4, [lowcut/nyq, highcut/nyq], btype='band')
        sub_data = sps.filtfilt(b, a, sub_data, axis=0)

        for i in range(sub_data.shape[1]):
            std_i = np.std(sub_data[:, i])
            if std_i > 0:
                sub_data[:, i] = sub_data[:, i] / std_i
        sub_data = np.clip(sub_data, -2.5, 2.5)

        fig, ax = plt.subplots(figsize=(14, 8))
        extent = [offsets_km[0], offsets_km[-1], plot_time[-1], plot_time[0]]
        im = ax.imshow(sub_data.T, aspect='auto', extent=extent,
                       cmap='seismic', interpolation='bilinear', vmin=-1.5, vmax=1.5)

        ref_station = station or (list(self.station_coordinates.keys())[0] if self.station_coordinates else None)
        if ref_station is not None:
            try:
                rec_lat = self.station_coordinates[ref_station]["latitude"]
                rec_lon = self.station_coordinates[ref_station]["longitude"]
                arrivals_p, arrivals_s = get_travel_times(event, self.get_current_params(), rec_lat, rec_lon)
                src = self.sources.get(event, {})
                origin_time = UTCDateTime(src.get("start_time"))
                p_rel = (origin_time + arrivals_p[0].time) - window_start_utc if arrivals_p else None
                s_rel = (origin_time + arrivals_s[0].time) - window_start_utc if arrivals_s else None
                epi_dist_km = degrees2kilometers(locations2degrees(
                    src.get("latitude"), src.get("longitude"), rec_lat, rec_lon))
                t_rel = (origin_time + epi_dist_km / 1.5) - window_start_utc

                for label, t_val, color in [("P", p_rel, 'red'), ("S", s_rel, 'blue'), ("T", t_rel, 'green')]:
                    if t_val is not None and plot_time.min() <= t_val <= plot_time.max():
                        ax.axhline(t_val, color=color, linestyle='--', linewidth=1.5, alpha=0.7)
                        ax.text(start_offset_km + 0.2, t_val, label, color=color,
                                fontsize=12, fontweight='bold', va='center', ha='left')
                        ax.text(end_offset_km - 0.2, t_val, label, color=color,
                                fontsize=12, fontweight='bold', va='center', ha='right')
            except Exception as e:
                logger.error(f"Unable to calculate theoretical phases: {e}")

        ax.set_xlabel("Offset [km]", fontsize=12, fontweight='bold')
        ax.set_ylabel("Time [s]", fontsize=12, fontweight='bold')
        ax.set_title(f"{event} - Band-passed 0.5-20 Hz", fontsize=14, pad=15)
        ax.grid(True, linestyle=':', alpha=0.3)
        cbar = plt.colorbar(im, ax=ax, label='Normalized Amplitude')
        cbar.ax.tick_params(labelsize=10)

        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        event_safe = event.replace(" ", "_")
        filename = os.path.join(save_path, f"{event_safe}{suffix}_record_section_2d.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        logger.info(f"2D plot saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_seismic_zoom(self, start_offset_km=30.8, end_offset_km=50.0,
                          t_min=5, t_max=30, event=None, suffix="", **kwargs):
        """
        Zoomed version of the seismic plot to view P and S arrivals.
        """

        if self.dfdas is None:
            logger.warning("No data in memory to plot.")
            return

        event = event or self.event
        logger.info(f"»»» Generating ZOOM Record Section for {event}...")

        data = self.dfdas.to_numpy()
        dt = self.dfdas.meta['dt']
        fs = 1.0 / dt
        time_axis = (self.dfdas.index - self.dfdas.index[0]).total_seconds()
        
        dx_native = self.dfdas.meta.get('dx', 1.0213)
        absolute_offsets_km = (self.dfdas.columns.to_numpy() * dx_native) / 1000.0

        spatial_mask = (absolute_offsets_km >= start_offset_km) & (absolute_offsets_km <= end_offset_km)
        sub_data = data[:, spatial_mask]
        offsets_km = absolute_offsets_km[spatial_mask]

        time_mask = (time_axis >= t_min) & (time_axis <= t_max)
        sub_data = sub_data[time_mask, :]
        plot_time = time_axis[time_mask]

        lowcut, highcut = 0.5, 20.0
        nyq = 0.5 * fs
        b, a = sps.butter(4, [lowcut/nyq, highcut/nyq], btype='band')
        sub_data = sps.filtfilt(b, a, sub_data, axis=0)

        for i in range(sub_data.shape[1]):
            std_i = np.std(sub_data[:, i])
            if std_i > 0:
                sub_data[:, i] = sub_data[:, i] / std_i
        sub_data = np.clip(sub_data, -2.5, 2.5)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        extent = [offsets_km[0], offsets_km[-1], plot_time[-1], plot_time[0]]
        im1 = ax1.imshow(sub_data.T, aspect='auto', extent=extent,
                         cmap='seismic', interpolation='bilinear', vmin=-1.5, vmax=1.5)
        ax1.set_xlabel("Offset [km]", fontsize=12)
        ax1.set_ylabel("Time [s]", fontsize=12)
        ax1.set_title(f"{event} - Zoom P/S arrivals", fontsize=14)
        ax1.grid(True, linestyle=':', alpha=0.3)

        mean_trace = np.mean(sub_data, axis=1)
        ax2.plot(mean_trace, plot_time, color='black', linewidth=1.5)
        ax2.set_xlabel("Mean Amplitude", fontsize=12)
        ax2.set_ylabel("Time [s]", fontsize=12)
        ax2.set_title("Stacked Trace", fontsize=14)
        ax2.grid(True, linestyle=':', alpha=0.3)
        
        # Linee teoriche
        window_start_utc = UTCDateTime(self.dfdas.meta['time'])
        ref_station = list(self.station_coordinates.keys())[0] if self.station_coordinates else None
        if ref_station is not None:
            try:
                rec_lat = self.station_coordinates[ref_station]["latitude"]
                rec_lon = self.station_coordinates[ref_station]["longitude"]
                src = self.sources.get(event, {})
                origin_time = UTCDateTime(src.get("start_time"))
                arrivals_p, arrivals_s = get_travel_times(event, self.get_current_params(), rec_lat, rec_lon)
                
                if arrivals_p:
                    p_rel = (origin_time + arrivals_p[0].time) - window_start_utc
                    ax1.axhline(y=p_rel, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.axhline(y=p_rel, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.text(0.1, p_rel, 'P', color='red', fontweight='bold')
                if arrivals_s:
                    s_rel = (origin_time + arrivals_s[0].time) - window_start_utc
                    ax1.axhline(y=s_rel, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.axhline(y=s_rel, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.text(0.1, s_rel, 'S', color='blue', fontweight='bold')
            except Exception as e:
                logger.error(f"Error computing phases: {e}")

        plt.tight_layout()
        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        filename = os.path.join(save_path, f"{event}_record_section_zoom.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        logger.info(f"Zoom plot saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_seismic_section(self, channel_step=1, scale_factor=1.0, start_offset_km=30.8,
                              end_offset_km=None, dx_meters=12.5, max_time_s=45.0,
                              station=None, event=None, suffix="", use_agc=True,
                              agc_window_s=2.0, agc_clip=4.0, **kwargs):
        """
        Record section with wiggle traces, band-passed 0.5-20 Hz, with P/S/T phases.
        """
        if self.dfdas is None:
            logger.warning("No data in memory to plot.")
            return

        event = event or self.event
        logger.info(f"»»» Generating Record Section for {event}...")

        data = self.dfdas.to_numpy()
        dt = self.dfdas.meta['dt']
        fs = 1.0 / dt
        window_start_utc = UTCDateTime(self.dfdas.meta['time'])
        time_axis = (self.dfdas.index - self.dfdas.index[0]).total_seconds()

        dx_native = self.dfdas.meta.get('dx', 1.0213001907746815)
        absolute_offsets_km = (self.dfdas.columns.to_numpy() * dx_native) / 1000.0

        if end_offset_km is None:
            end_offset_km = start_offset_km + 21.5

        logger.info(f"-> Diagnostics: dx={dx_native:.4f} m/channel | dt={dt:.5f}s (fs={fs:.1f}Hz) | "
              f"available channel coverage: {absolute_offsets_km.min():.2f} - {absolute_offsets_km.max():.2f} km | "
              f"data window: {time_axis.min():.1f} - {time_axis.max():.1f} s | "
              f"window start (UTC): {window_start_utc}")

        spatial_mask = (absolute_offsets_km >= start_offset_km) & (absolute_offsets_km <= end_offset_km)
        sub_data = data[:, spatial_mask]
        offsets_km = absolute_offsets_km[spatial_mask]

        if sub_data.shape[1] == 0:
            logger.error("ERROR: no channels in requested offset range.")
            return

        time_mask = time_axis <= max_time_s
        sub_data = sub_data[time_mask, :]
        plot_time = time_axis[time_mask]

        lowcut, highcut = 0.5, 20.0
        nyq = 0.5 * fs
        b, a = sps.butter(4, [lowcut/nyq, highcut/nyq], btype='band')
        sub_data = sps.filtfilt(b, a, sub_data, axis=0)

        trace_spacing_km = dx_meters / 1000.0
        target_offsets = np.arange(start_offset_km, end_offset_km, trace_spacing_km)
        final_offsets = []
        final_data_indices = []
        for target in target_offsets:
            idx = np.abs(offsets_km - target).argmin()
            if idx not in final_data_indices:
                final_data_indices.append(idx)
                final_offsets.append(offsets_km[idx])

        final_data = sub_data[:, final_data_indices]
        final_offsets = np.array(final_offsets)
        logger.info(f"-> Plotting {final_data.shape[1]} traces (target spacing: {dx_meters} m).")

        fig, ax = plt.subplots(figsize=(15, 8.5))
        max_deflection_km = trace_spacing_km * 0.45 * scale_factor
        agc_win_samples = max(1, int(agc_window_s * fs))
        agc_kernel = np.ones(agc_win_samples) / agc_win_samples

        # One LineCollection for all traces (instead of one ax.plot per channel).
        # Shared time axis goes in once; x is filled per trace below.
        segments = np.empty((final_data.shape[1], plot_time.size, 2))
        segments[:, :, 1] = plot_time
        n_seg = 0

        for i in range(final_data.shape[1]):
            traccia = final_data[:, i]
            traccia = sps.detrend(traccia)
            traccia = traccia - np.median(traccia)
            edge_samples = int(0.5 * fs)
            if edge_samples > 0:
                traccia[:edge_samples] = 0
                traccia[-edge_samples:] = 0

            if use_agc:
                local_power = sps.fftconvolve(traccia**2, agc_kernel, mode='same')
                local_rms = np.sqrt(np.maximum(local_power, 0))
                floor = np.median(local_rms[local_rms > 0]) * 0.05 if np.any(local_rms > 0) else 1e-12
                local_rms = np.maximum(local_rms, floor if floor > 0 else 1e-12)
                traccia_norm = np.clip(traccia / local_rms, -agc_clip, agc_clip)
                traccia_km = (traccia_norm / agc_clip) * max_deflection_km
            else:
                std_local = np.std(traccia)
                if std_local == 0:
                    continue
                traccia_norm = np.clip(traccia / std_local, -2.5, 2.5)
                traccia_km = (traccia_norm / 2.5) * max_deflection_km

            segments[n_seg, :, 0] = final_offsets[i] + traccia_km
            n_seg += 1

        if n_seg:
            ax.add_collection(LineCollection(segments[:n_seg], colors='#2c3e50',
                                             linewidths=0.4, alpha=0.8))

        # Fasi teoriche
        ref_station = station or (list(self.station_coordinates.keys())[0] if self.station_coordinates else None)
        if ref_station is not None:
            try:
                rec_lat = self.station_coordinates[ref_station]["latitude"]
                rec_lon = self.station_coordinates[ref_station]["longitude"]
                arrivals_p, arrivals_s = get_travel_times(event, self.get_current_params(), rec_lat, rec_lon)
                src = self.sources.get(event, {})
                origin_time = UTCDateTime(src.get("start_time"))
                p_rel = (origin_time + arrivals_p[0].time) - window_start_utc if arrivals_p else None
                s_rel = (origin_time + arrivals_s[0].time) - window_start_utc if arrivals_s else None
                epi_dist_km = degrees2kilometers(locations2degrees(
                    src.get("latitude"), src.get("longitude"), rec_lat, rec_lon))
                t_rel = (origin_time + epi_dist_km / 1.5) - window_start_utc

                for label, t_val in [("P", p_rel), ("S", s_rel), ("T", t_rel)]:
                    if t_val is not None and plot_time.min() <= t_val <= plot_time.max():
                        ax.axhline(t_val, color='crimson', linestyle='--', linewidth=0.7, alpha=0.7)
                        ax.text(start_offset_km - 0.3, t_val, label, color='crimson',
                                fontsize=10, fontweight='bold', va='center', ha='right')
                        ax.text(end_offset_km + 0.3, t_val, label, color='crimson',
                                fontsize=10, fontweight='bold', va='center', ha='left')
                    else:
                        logger.warning(f"WARNING: phase {label} at t={t_val} is outside the plotted window "
                              f"({plot_time.min():.1f}-{plot_time.max():.1f}s).")
            except Exception as e:
                logger.error(f"Unable to calculate theoretical phases: {e}")

        ax.set_xlim(start_offset_km, end_offset_km)
        ax.set_ylim(max_time_s, 0.0)
        ax.set_xlabel("Offset [km]", fontsize=10, labelpad=10)
        ax.set_ylabel("Time [s]", fontsize=10, labelpad=10)
        ax.set_title(f"{event} record section", fontsize=12, pad=15)
        ax.grid(True, linestyle=':', color='gainsboro', linewidth=0.5)
        ax.tick_params(top=True, bottom=True, left=True, right=True,
                       labeltop=True, labelbottom=True, labelleft=True, labelright=True, labelsize=9)
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        event_safe = event.replace(" ", "_")
        filename = os.path.join(save_path, f"{event_safe}_record_section.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        logger.info(f"Figure saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_seismic_enhanced(self, channel_step=8, scale_factor=15.0, 
                               start_offset_km=30.8, end_offset_km=50.0,
                               dx_meters=102.0, max_time_s=40.0,
                               station=None, event=None, suffix="", **kwargs):
        """
        Enhanced version of the seismic plot optimized for the M2.7 Desertas earthquake,
        as in Figure 7 of the paper.
        """

        if self.dfdas is None:
            logger.warning("No data in memory to plot.")
            return

        event = event or self.event
        logger.info(f"»»» Generating ENHANCED Record Section for {event}...")

        data = self.dfdas.to_numpy()
        dt = self.dfdas.meta['dt']
        fs = 1.0 / dt
        window_start_utc = UTCDateTime(self.dfdas.meta['time'])
        time_axis = (self.dfdas.index - self.dfdas.index[0]).total_seconds()

        dx_native = self.dfdas.meta.get('dx', 1.0213)
        absolute_offsets_km = (self.dfdas.columns.to_numpy() * dx_native) / 1000.0

        spatial_mask = (absolute_offsets_km >= start_offset_km) & (absolute_offsets_km <= end_offset_km)
        sub_data = data[:, spatial_mask]
        offsets_km = absolute_offsets_km[spatial_mask]

        if sub_data.shape[1] == 0:
            logger.error("ERROR: no channels in the requested range.")
            return

        time_mask = time_axis <= max_time_s
        sub_data = sub_data[time_mask, :]
        plot_time = time_axis[time_mask]

        # 0.5-20 Hz filter
        lowcut, highcut = 0.5, 20.0
        nyq = 0.5 * fs
        b, a = sps.butter(4, [lowcut/nyq, highcut/nyq], btype='band')
        sub_data = sps.filtfilt(b, a, sub_data, axis=0)

        # Spatial decimation
        step = max(1, channel_step)
        sub_data = sub_data[:, ::step]
        offsets_km = offsets_km[::step]

        # AGC normalization
        for i in range(sub_data.shape[1]):
            window_samples = int(2.0 * fs)
            kernel = np.ones(window_samples) / window_samples
            power = sps.fftconvolve(sub_data[:, i]**2, kernel, mode='same')
            rms = np.sqrt(np.maximum(power, 1e-12))
            rms[rms < np.median(rms) * 0.1] = np.median(rms) * 0.1
            sub_data[:, i] = np.clip(sub_data[:, i] / rms, -3, 3)

        # Create figure with 3 panels
        fig, axes = plt.subplots(3, 1, figsize=(14, 12))

        # 1. Wiggle plot
        ax1 = axes[0]
        trace_spacing_km = offsets_km[1] - offsets_km[0] if len(offsets_km) > 1 else 0.1
        max_deflection = trace_spacing_km * 0.6 * scale_factor

        # One LineCollection for all traces (instead of one ax.plot per channel)
        segments = np.empty((sub_data.shape[1], plot_time.size, 2))
        segments[:, :, 1] = plot_time
        for i in range(sub_data.shape[1]):
            trace = sub_data[:, i]
            trace = sps.detrend(trace)
            trace = trace - np.median(trace)
            segments[i, :, 0] = offsets_km[i] + trace * max_deflection / 3.0
        ax1.add_collection(LineCollection(segments, colors='k',
                                          linewidths=0.3, alpha=0.7))

        ax1.set_ylabel("Time [s]")
        ax1.set_xlim(start_offset_km, end_offset_km)
        ax1.set_ylim(max_time_s, 0)
        ax1.set_title("Wiggle Plot (AGC normalized)")
        ax1.grid(True, alpha=0.2)

        # 2. Heatmap 2D
        ax2 = axes[1]
        extent = [offsets_km[0], offsets_km[-1], plot_time[-1], plot_time[0]]
        im = ax2.imshow(sub_data.T, aspect='auto', extent=extent,
                        cmap='seismic', interpolation='bilinear',
                        vmin=-1.5, vmax=1.5)
        ax2.set_ylabel("Time [s]")
        ax2.set_xlabel("Offset [km]")
        ax2.set_title("Heatmap (AGC normalized)")
        plt.colorbar(im, ax=ax2, label='Normalized Amplitude')

        # 3. Stacked trace with phases
        ax3 = axes[2]
        stack = np.mean(sub_data, axis=1)
        ax3.plot(stack, plot_time, 'k-', linewidth=1.5)
        ax3.set_xlabel("Stacked Amplitude")
        ax3.set_ylabel("Time [s]")
        ax3.set_title("Stacked Trace")
        ax3.grid(True, alpha=0.3)
        ax3.set_ylim(max_time_s, 0)

        # Aggiungi fasi teoriche
        ref_station = station or (list(self.station_coordinates.keys())[0] if self.station_coordinates else None)
        if ref_station is not None:
            try:
                rec_lat = self.station_coordinates[ref_station]["latitude"]
                rec_lon = self.station_coordinates[ref_station]["longitude"]
                src = self.sources.get(event, {})
                origin_time = UTCDateTime(src.get("start_time"))
                arrivals_p, arrivals_s = get_travel_times(event, self.get_current_params(), rec_lat, rec_lon)
                
                if arrivals_p:
                    p_abs = origin_time + arrivals_p[0].time
                    p_rel = p_abs - window_start_utc
                    ax1.axhline(y=p_rel, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.axhline(y=p_rel, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax3.axhline(y=p_rel, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax3.text(0.1, p_rel, 'P', color='red', fontweight='bold')
                if arrivals_s:
                    s_abs = origin_time + arrivals_s[0].time
                    s_rel = s_abs - window_start_utc
                    ax1.axhline(y=s_rel, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax2.axhline(y=s_rel, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax3.axhline(y=s_rel, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
                    ax3.text(0.1, s_rel, 'S', color='blue', fontweight='bold')
                # T-phase
                epi_dist_km = degrees2kilometers(locations2degrees(
                    src.get("latitude"), src.get("longitude"), rec_lat, rec_lon))
                t_abs = origin_time + epi_dist_km / 1.5
                t_rel = t_abs - window_start_utc
                if 0 <= t_rel <= max_time_s:
                    ax1.axhline(y=t_rel, color='green', linestyle='--', linewidth=1.0, alpha=0.5)
                    ax2.axhline(y=t_rel, color='green', linestyle='--', linewidth=1.0, alpha=0.5)
                    ax3.axhline(y=t_rel, color='green', linestyle='--', linewidth=1.0, alpha=0.5)
                    ax3.text(0.1, t_rel, 'T', color='green', fontweight='bold')
            except Exception as e:
                logger.error(f"Error computing phases: {e}")

        plt.suptitle(f"{event} - Seismic Record Section (M2.7 Desertas)", fontsize=14, fontweight='bold')
        plt.tight_layout()

        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        filename = os.path.join(save_path, f"{event}_enhanced_section.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        logger.info(f"Enhanced plot saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_record_sections(self, dataset=None, event=None, suffix=""):
        # Use class attributes for configuration to integrate with DASProcessor
        ep, pp, fp = self.exec_params, self.path_params, self.freq_params
        station_coordinates = self.station_coordinates

        event_info = self.sources.get(event, {})
        station=ep["station"]
    
        logger.info("Plotting record section.")
        result=subprocess.run(["bash","./plot_GMTmap.sh",
                                      str(event_info.get("longitude")),
                                      str(event_info.get("latitude")),
                                      f'"{str(event_info.get("magnitude"))} {event_info.get("depth")}km"',
                                      f'{str(event)}_map.ps'],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
        filename=f'{pp.get("figures_path")}/{event}_record-section{f"_{suffix}" if suffix else ""}'
    
        dist_az_baz=gps2dist_azimuth(event_info['latitude'],event_info['longitude'],
                                     self.station_coordinates.get(station)["latitude"],
                                     self.station_coordinates.get(station)["longitude"])
        if dist_az_baz[0]<100000:
            local = True
            logger.info("Distance to event is less than 100 km. Considering it a local event for plotting.")
        else:
            local = False
            logger.info("Distance to event is greater than 100 km. Considering it a regional/teleseismic event.")

        if ep.get("verbose"):
            logger.debug("Standard output:\n%s", result.stdout)
            logger.debug("Standard error:\n%s", result.stderr)
            if result.returncode == 0:
                logger.debug("Bash script executed successfully.")
            else:
                logger.warning("Bash script returned an error (return code %d).", result.returncode)
    
        plt.rcParams['figure.dpi'] = 1200
        fig=plt.figure(clear=True,figsize=(10,6))
    
    
    # get land stations
        
        source_latitude = event_info.get("latitude")
        source_longitude = event_info.get("longitude")
        source_depth = event_info.get("depth")
        station=ep.get("land_station")
        station_lat=station_coordinates.get(station,{}).get("latitude")
        station_lon=station_coordinates.get(station,{}).get("longitude")
    
        i91 = TauPyModel(model="iasp91")
    
        #filtering according to distance
        if local:
            max_freq=30
            min_freq=.5
        else:
            max_freq=.1
            min_freq=0.01
    
        fig = plt.figure(figsize=(40, 8))  #12,4 #54,8  #12,4
    
        if ep.get("plot_stack"):
            fig = plt.figure(figsize=(40, 8))  #12,4 #54,8  #12,4
            gs = GridSpec(1, 5, width_ratios=[9, .7, .7, .7, .7])
        elif ep.get("plot_land_station"):
            fig = plt.figure(figsize=(40, 8))  #12,4 #54,8  #12,4
            gs = GridSpec(1, 4, width_ratios=[9, .7, .7, .7])
        else:
            fig = plt.figure(figsize=(30, 8))  #12,4 #54,8  #12,4
            gs = GridSpec(1, 1)
    
        arrivals_p = i91.get_travel_times_geo(source_depth_in_km=source_depth, receiver_latitude_in_deg=station_lat,\
                                                  receiver_longitude_in_deg=station_lon, source_latitude_in_deg=source_latitude,\
                                                  source_longitude_in_deg=source_longitude,phase_list=('p','sP','s','P','PP','S','SS','SKS'), ray_param_tol=0.1)
        
        # calculate trim points    
        time_delta=event_info["window_in"]+(UTCDateTime(event_info["end_time"])-UTCDateTime(event_info["start_time"]))
        
        pdiff_time=UTCDateTime(event_info.get("start_time"))+arrivals_p[0].time
        start_trim_time=pdiff_time-datetime.timedelta(seconds=event_info.get("window_in"))
        end_trim_time=start_trim_time+datetime.timedelta(seconds=time_delta)
    
        date_format = DateFormatter("%d-%m-%Y %H:%M:%S")
    
        # Reference offset relative to event origin time
        p_time = arrivals_p[0].time - event_info.get('window_in')
        event_time = UTCDateTime(event_info.get('start_time'))
    
        for arr in arrivals_p:
            # Keep time as float offset (seconds)
            arr.time = round(arr.time, 3)
            # Print absolute arrival time
            logger.debug("%s %s", arr.name, event_time + arr.time)
            # Normalise relative arrival time
            arr.time -= p_time
            # Print updated absolute arrival time
            logger.debug("%s %s", arr.name, event_time + arr.time)
    
        # get DAS data

        # If no dataset passed, try to use the class-loaded self.dfdas via the class method.
        # If a non-Stream `dataset` (e.g., pandas DataFrame) is provided, temporarily set
        # `self.dfdas` so we can reuse `self.dfdas_to_mseed` and avoid importing the
        # standalone helper function here.
        restore_dfdas = False
        if dataset is None:
            # Prefer using already-loaded data inside the class
            st = None
            if hasattr(self, 'mseed_stream') and self.mseed_stream is not None:
                st = self.mseed_stream
            elif self.dfdas is not None:
                st = self.dfdas_to_mseed(event, write=False)
        elif isinstance(dataset, Stream):
            st = dataset
        else:
            # dataset is likely a pandas DataFrame (dfdas). Temporarily assign and convert.
            orig_dfdas = getattr(self, 'dfdas', None)
            self.dfdas = dataset
            restore_dfdas = True
            st = self.dfdas_to_mseed(event, write=False)

        if restore_dfdas:
            self.dfdas = orig_dfdas

        # If conversion didn't produce a Stream, bail out
        if st is None:
            logger.warning("No Stream available for plotting.")
            return

        st.filter("bandpass", freqmin=min_freq, freqmax=max_freq, corners=4, zerophase=True)
        
        # decimate DAS for plotting
        # Loop through each trace in the Stream
        # Create a new Stream object containing only every other trace
        every_other_trace = []
    
        if local:
            decim=400 #2000 #4000 #2000
        else:
            decim=100
        if len(st)> decim:
            for i, trace in enumerate(st):
            # Check if the index is even (i.e., every other trace)
                if i % int(len(st)/decim) == 0:
                    every_other_trace.append(trace)
            st_every_other = Stream(traces=every_other_trace)
        else:
            st_every_other = st.copy()            
    
        st_event=Stream()
        if ep.get("plot_land_station"):
            logger.info(f'Getting channels {station_coordinates.get(station,{}).get("channels")} from {ep.get("land_station")}')
            if "MSEED_file" in self.sources.get(event,{}):
                st_land=read(os.path.join(pp.get("MSEED_path", ""),self.sources.get(event, {}).get("MSEED_file", "")))
            else:
                st_land=read(pp.get("mseed_land_stations"))
            st_event=st_land.select(channel=station_coordinates.get(station,{}).get("channels"),station=ep.get("land_station"))
            st_event.trim(starttime=start_trim_time, endtime=end_trim_time)
        
        if len(st_event)>0:
            if start_trim_time<st_event[0].stats.starttime-0.01:
                logger.warning(f'{st_event[0].stats.station}.{st_event[0].stats.channel} has not enough data before the event.')
                logger.info(f'{st_event[0].stats.starttime}>{start_trim_time} diff={start_trim_time-st_event[0].stats.starttime}')
                # window_correction=start_trim_time-st_event[0].stats.starttime
                # sources.get(event, {})["window_in"] = sources.get(event, {}).get("window_in")+window_correction
                # re_trim=True
            if end_trim_time<st_event[0].stats.endtime-0.01:
                logger.warning(f'{st_event[0].stats.station}.{st_event[0].stats.channel} has not enough data after the event.')
                logger.info(f'{st_event[0].stats.endtime}>{end_trim_time} diff={end_trim_time-st_event[0].stats.endtime}')
                # window_correction=end_trim_time-st_event[0].stats.endtime
                # sources.get(event, {})["window_out"] = sources.get(event, {}).get("window_out")+window_correction
            # re_trim=True
        else:
            logger.warning(f'The selected channel ({station_coordinates.get(station,{}).get("channels")}) doesn\'t exist in {os.path.join(pp.get("MSEED_path", ""),self.sources.get(event, {}).get("MSEED_file", ""))}.')
            # Try to download land station data from a configured EIDA node if available
            try:
                logger.info("Attempting to download land station data from configured EIDA node...")
                st_dl = self.download_waveform_data(station_code=station, starttime=start_trim_time, endtime=end_trim_time)
                if st_dl is not None and len(st_dl) > 0:
                    # Select requested channels and trim to window
                    try:
                        st_event = st_dl.select(channel=station_coordinates.get(station, {}).get("channels"), station=station)
                        st_event.trim(starttime=start_trim_time, endtime=end_trim_time)
                        logger.info(f"Successfully downloaded and trimmed data for station {station}.")
                    except Exception as e:
                        logger.warning(f"Downloaded data present but required channels not found: {e}")
                else:
                    logger.warning("No data returned from EIDA for the requested interval/station.")
            except Exception as e:
                logger.error(f"EIDA download attempt failed: {e}")
    
        if len(st_event)>0:
            unique_phases=[]
            st_event.detrend()
            st_event.taper(max_percentage=0.05, type='hann')
            st_event.filter("bandpass", freqmin=fp.get("land_station_min_freq"), freqmax=fp.get("land_station_max_freq"), corners=8, zerophase=True)
            logger.info(f'Land stations bandpass filtered between {fp.get("land_station_min_freq")} and {fp.get("land_station_max_freq")} Hz')
            ylim=(round(max(abs(st_event[0].data))/5)+1)*5
    
            if len(st_event)>3:
                st_z = st_event.select(channel='*Z')[:1]
                st_n = st_event.select(channel='*N')[:1]
                st_e = st_event.select(channel='*E')[:1]
                st_event.clear()
    
                # Append the selected traces back into st_event
                st_event += st_z
                st_event += st_n
                st_event += st_e
    
        if ep.get("plot_stack"):
            trace_ranges={(0,10),(20,30),(50,60),(- 10,None)}
            for start,end in trace_ranges:
                last_traces=st_every_other[start:end]
                st_stack=Stream()
                st_stack.extend(last_traces)
                st_event +=st_stack.stack()
                st_event[-1].stats.station="DAS"
                st_event[-1].stats.location=f'{last_traces[0].stats.distance}-{last_traces[-1].stats.distance}'
            ylim_stack=1.2e-8
    
        if len(st_event)>0:
            ax1 = [None] * len(st_event)
            ax2 = [None] * len(st_event)
            for idx,event_channel in enumerate(st_event): #st_event
                ax1[idx] = fig.add_subplot(gs[idx+1])
                ax1[idx].set_title(f'{event_channel.stats.station}-{event_channel.stats.channel}')
                filename+=f'{st_event[idx].stats.station}-{st_event[idx].stats.channel}'
                
                ax2[idx] = ax1[idx].twinx()
                ax2[idx].yaxis.set_major_formatter(date_format)
                ax2[idx].plot(event_channel.data, event_channel.times("matplotlib"), color="black", linewidth=0.05)
                ax1[idx].plot([0, 0], [0, time_delta], color="white")
                ax2[idx].invert_yaxis()
                ax1[idx].margins(0, 0)
                ax2[idx].margins(0, 0)
                ax1[idx].invert_yaxis()
                logger.info(f'Adding {event_channel.stats.station}-{event_channel.stats.channel}')
                if event_channel.stats.station=="DAS":
                    ax1[idx].set_xlim(-ylim_stack, ylim_stack)
                    ax1[idx].set_title(f'{event_channel.stats.station}-{event_channel.stats.location}')
                else:
                    ax1[idx].set_xlim(-ylim, ylim)
                    ax1[idx].set_title(f'{event_channel.stats.station}-{event_channel.stats.channel}')
                ax1[idx].set_ylim(time_delta,0)
            
            label_space = 1 if local else 10
            prev=0
            unique_phases=[]
            if not ep.get("plot_land_station"):
                ylim=ylim_stack
            for arr in arrivals_p:
                if arr.name not in unique_phases:
                    unique_phases.append(arr.name)
                    if arr.time > prev+label_space and arr.time < time_delta:
                        prev=arr.time
                        plot_time=mdates.date2num(UTCDateTime(self.sources.get(event, {}).get("time"))+arr.time)
                        logger.debug("%s %s %s %s %s", "UNIQUE ARRIVAL", arr.name, arr.time, plot_time, UTCDateTime(self.sources.get(event, {}).get("time"))+arr.time)
                        ax2[0].axhline(y=plot_time, color='red')
                        ax2[0].annotate(arr.name, xy=(ylim*0.9, plot_time), xytext=(-100, 0),
                                    textcoords='offset pixels')
            logger.debug("%s %s", "Unique arrivals", unique_phases)
            
        ax0 = fig.add_subplot(gs[0])
    
        # calculate arrival times for every trace of the cable
        offset_traveltime_arrival = {}
        for tr in st_every_other:
            part_arrival=[]
            station_lat=tr.stats.coordinates.latitude
            station_lon=tr.stats.coordinates.longitude
            arrivals = i91.get_travel_times_geo(source_depth_in_km=source_depth,
                                                receiver_latitude_in_deg=station_lat,
                                                receiver_longitude_in_deg=station_lon,
                                                source_latitude_in_deg=source_latitude,
                                                source_longitude_in_deg=source_longitude,
                                                phase_list=('p','sP','s','P','PP','S','SS','SKS'),
                                                ray_param_tol=0.1)

            logger.debug("%s %s %s", event_info["start_time"], arrivals[0], self.clock_offset)
            pdiff_time=UTCDateTime(event_info.get("start_time"))+arrivals[0].time+datetime.timedelta(seconds=self.clock_offset)
            start_trim_time=pdiff_time-datetime.timedelta(seconds=event_info.get("window_in"))
        
            p_time=arrivals[0].time-event_info.get("window_in")
        
            unique_phases=[]
            for arrival in arrivals:
                if arrival.name not in unique_phases:
                    unique_phases.append(arrival.name)
                    arrival.time=round(arrival.time, 3)
                    arrival.time-=p_time-1
                    distance_km = tr.stats.distance / 1000
                    if arrival.name not in offset_traveltime_arrival:
                        offset_traveltime_arrival[arrival.name] = []
                    offset_traveltime_arrival[arrival.name].append([distance_km,arrival.time])
    
        if ep.get("plot_section_colour"):
            st_every_other.plot(fig=fig, type='section', norm_method="trace", time_down=True, orientation='vertical',linewidth=0.1,fillcolors=('blue','red'), scale=2.5, dpi=900)
        else:
            st_every_other.plot(fig=fig, type='section', norm_method="trace", time_down=True, orientation='vertical',linewidth=0.1, scale=1.5, dpi=900)
    
        ax0.set_title(f'{event} M{str(self.sources.get(event, {}).get("magnitude"))}')
        ax3=ax0.twinx()
        #ax3=ax0.twiny()
        ax3.sharey(ax0)
        if len(st_event)>0:
            ax1[0].sharey(ax0)
            yticklabels = ax1[0].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        # ax3.invert_yaxis()
        # ax3.set_ylim(time_delta-0.0001,0)
        prev=0
        annot=""
        label_space = 1 if local else 10
        logger.debug("%s %s %s", "LABEL SPACE:", local, label_space)
    
        for name, data in offset_traveltime_arrival.items():
            x=[point[0] for point in data]
            y=[point[1] for point in data]
            ax3.plot(x,y, linestyle='-',linewidth=.5, color='red')
            annot=annot+" "+name
            first_point=data[0]
            if first_point[1] > prev+label_space and first_point[1] < time_delta:
                ax3.annotate(annot,xy=(first_point[0],first_point[1]),textcoords="offset points",xytext=(-5,0),ha='center',fontsize=8)
                annot=""
            prev=first_point[1]
        
        plt.tight_layout()
        plt.savefig(filename+'.eps',dpi=900, format='eps')
        plt.savefig(filename+'.pdf',dpi=300, format='pdf')
        #plt.savefig(filename+'.png',dpi=900, format='png')
        # Add the existing figure as a subplot
        plt.show()
        plt.close()
    
        logger.debug(f'Interspace spacing={st_every_other[1].stats.distance-st_every_other[0].stats.distance} from {st_every_other[0].stats.distance} to {st_every_other[-1].stats.distance}. {len(st_every_other)} traces.')
        logger.info("Finished plotting record section.")
