import inspect
import datetime,time
import os
import numpy as np
import warnings
from collections.abc import Iterable
from sympy import symbols, sympify
import re
import xdas as xd
import builtins
import string
import math
import h5py
import pandas as pd
import pickle
import json
import logging
from pathlib import Path

from contextlib import contextmanager

from obspy.core.util import AttribDict
from obspy.taup import TauPyModel
from obspy import UTCDateTime, read, Trace, Stream
from obspy.io.segy.segy import SUFile,SEGYTraceHeader, SEGYBinaryFileHeader
from obspy.io.segy.core import _read_segy
from obspy.signal.array_analysis import array_processing,array_transff_wavenumber,array_transff_freqslowness,get_geometry
from obspy.signal.invsim import cosine_taper
from obspy.imaging.cm import obspy_sequential
from obspy.geodetics import degrees2kilometers,locations2degrees
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Comment, Response, Equipment, Operator
from obspy.core.inventory import InstrumentSensitivity, Person, ResponseStage, PolesZerosResponseStage, CoefficientsTypeResponseStage
from obspy.core.inventory import read_inventory
#from obspy.signal.invsim import corn_freq_2_paz

from matplotlib import mlab

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap # Normalize
from matplotlib.patches import Patch #,ConnectionPatch, Normalize, Rectangle,Patch
from matplotlib.transforms import (Bbox, TransformedBbox,
                                   blended_transform_factory)
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import (BboxConnector,
                                                   BboxConnectorPatch,
                                                   BboxPatch)
import scipy.signal as sps
from scipy import ndimage
from scipy import signal
from scipy.signal import hilbert, square, ShortTimeFFT
from scipy.signal.windows import gaussian,hann

import concurrent.futures
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
import subprocess
from multiprocessing import Process

from memory_profiler import profile

#check if code is running in jupyter
if hasattr(builtins, "__IPYTHON__"):
    from tqdm.notebook import tqdm
    from ipywidgets import FloatProgress
    #%matplotlib widget
    #%config InlineBackend.print_figure_kwargs = {'bbox_inches':None}
else:
    from tqdm import tqdm

#from defaults import offset_ranges,channel_ranges
from sources import sources

# define class for storing results

class DASProcessor:
    def __init__(self, global_params):
        self.exec_params = global_params.get('exec_params', {})
        self.path_params = global_params.get('path_params', {})
        self.sources = global_params.get('sources', {})
        self.channel_ranges = global_params.get('channel_ranges', {})
        self.offset_ranges = global_params.get('offset_ranges', {})
        self.station_coordinates = global_params.get('station_coordinates', {})
        self.acquisition_params = global_params.get('acquisition_params', {})
        self.trace_header = global_params.get('trace_header', {})
        self.freq_params = global_params.get('freq_params', {})

        self.stationxml_params={
                                "input_sampling_rate": [],
                                "filter_a_coefficients": [],
                                "filter_b_coefficients": [],
                                }
        self.network = self.exec_params.get("network", "5V")  # Default to '5V' if not specified

        # Internal state
        self.event_info = {}
        self.files = []
        self.dfdas = None
        self.mseed_stream = None
        self.filemeta = None
        self.roi_step = 0
        self.first_chan = 0
        self.last_chan = 0
        self.event= None
        self.filename_suffix = ""
        
        # New attributes to store parameters calculated during processing
        self.datetime_start = None
        self.duration = None
        self.channels = []
        self.figsize = (10, 8) # Default

    # Channel Results
        self.channels = np.array([])
        self.min_channel = []
        self.max_channel = []
        self.channel_step = []
        self.chIndex = []

        # Set up basic logging configuration with execution timestamp
        run_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = f'DAS_processor_{run_time}.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            filename=log_filename
        )

    def get_current_params(self):
        """Returns the internal state in the old dictionary format for compatibility."""
        return {
            'exec_params': self.exec_params,
            'path_params': self.path_params,
            'acquisition_params': self.acquisition_params,
            'sources': self.sources,
            'station_coordinates': self.station_coordinates,
            'freq_params': self.freq_params,
            'netwrk': self.network,
            'trace_header': self.trace_header[self.network],
            'stationxml_params': self.stationxml_params,
        }
    
    @contextmanager
    def temporary_event_paths(self, event_name):
        """
        Temporarily overrides paths and execution parameters for a specific event 
        and restores them automatically.
        """
        # 1. Store original values safely from their respective dictionaries
        original_paths = {k: self.path_params.get(k) for k in ["experiment_path", "MSEED_path", "events_path"]}
        original_exec = {"save_DFDAS": self.exec_params.get("save_DFDAS")}


        path_overrides = {
            "experiment_path": self.event_info.get("DAS_data_path"),
            "MSEED_path": self.event_info.get("mseed_output_path"),
            "events_path": self.event_info.get("events_path"),
        }
        exec_overrides = {
            "save_DFDAS": self.event_info.get("save_DAS_event")
        }

        try:
            # 3. Apply path overrides (checking explicitly for None so empty strings/False aren't skipped)
            for key, new_val in path_overrides.items():
                if new_val is not None:
                    print(f"Override: Setting path_params[{key}] to {new_val}")
                    self.path_params[key] = new_val

            # Apply execution overrides (using 'is not None' so False is successfully set!)
            for key, new_val in exec_overrides.items():
                if new_val is not None:
                    print(f"Override: Setting exec_params[{key}] to {new_val}")
                    self.exec_params[key] = new_val

            yield  # <--- Processing runs here

        finally:
            # 4. RESET to original defaults securely
            for key, old_val in original_paths.items():
                self.path_params[key] = old_val

            for key, old_val in original_exec.items():
                self.exec_params[key] = old_val

            print("Paths and execution parameters restored to defaults.")

    def run_process_chain(self, suffix=""):
        """
            Iterates through the 'process_chain' defined in the source
            and executes methods dynamically.
            Passes a suffix (e.g., '_timelapse_1') to methods for unique filenames.
            The stored chain definition is never modified, so this method can be
            called repeatedly with different suffixes.
        """

        chain = self.event_info.get('process_chain', [])
        ep = self.exec_params

        for step in chain:
            operation = step.get('operation')

            # Work on a copy: the original step['params'] stays untouched
            params = dict(step.get('params') or {})

            logging.info(f"Processing step: {operation} with params: {params}")
            
            if not hasattr(self, operation):
                print(f"Warning: Method {operation} not found in DASProcessor.")
                logging.warning(f"Method {operation} not found in DASProcessor.")
                continue

            method = getattr(self, operation)
            logging.info(f"Executing step: {operation} with params: {params}")
            if ep.get("verbose"):
                print(f"--- Executing: {operation} with params {params} ---")
            
            # 1. Accepted arguments of the method
            sig = inspect.signature(method)
            accepted_keys = sig.parameters.keys()
            accepts_kwargs = any(
                p.kind is inspect.Parameter.VAR_KEYWORD
                for p in sig.parameters.values()
            )

           # 2. Inject runtime values (only into the copy)
            if 'event' in accepted_keys and 'event' not in params:
                params['event'] = self.event

            # 'suffix' always comes from the current call, so a stale value
            # can never leak in from a previous run
            if 'suffix' in accepted_keys:
                params['suffix'] = suffix

            # 3. Filter to what the method accepts (unless it takes **kwargs)
            if accepts_kwargs:
                call_params = params
            else:
                call_params = {k: v for k, v in params.items() if k in accepted_keys}

            # 4. Call
            method(**call_params)
            logging.info(f"Completed step: {operation}")

    def get_time_lapse_windows(self, event):
        """Generates start times for continuous windows."""
        start_time = UTCDateTime(self.event_info.get("start_time"))
        end_time = UTCDateTime(self.event_info.get("end_time"))
        duration = self.event_info.get("window_duration")

        start_times = []
        current = start_time
        while current < end_time:
            start_times.append(current)
            current += duration

        print(f'{len(start_times)} windows of {duration}s to process')
        return start_times, duration
    
    def process_time_lapse(self, event):
        """Main loop for time-lapse processing."""
        print(f'»»» Working on time lapse: {event}')
        
        # 1. Generate the windows
        start_times, window_duration = self.get_time_lapse_windows(event)
        self.duration = window_duration # Set window length for the processor

        events_path = self.path_params.get("events_path")
        event_safe = event.replace(" ", "_")

        # 2. Iterate through windows
        for i, start_dt in enumerate(tqdm(start_times, desc="Time lapse windows")):
            suffix = f"_{start_dt.strftime('%Y%m%d_%H%M')}"
            target_hdf5 = os.path.join(events_path, f"{event_safe}{suffix}.hdf5")

            # Check if we should skip before loading raw data
            if os.path.exists(target_hdf5):
                print(f"Window {i+1}/{len(start_times)}: Skipping (Already exists: {os.path.basename(target_hdf5)})")
                continue

            print(f"\n--- Processing Window {i+1}/{len(start_times)}: {start_dt} ---")
            
            # Update the start time for this iteration
            self.datetime_start = start_dt
            
            # 3. Load data for this specific window
            # This triggers check_data_access -> get_roi_step -> actual loading
            success = self.load_event_data(event, mode="time_lapse")
            
            if success:
                # 4. Execute the chain for this window
                # We pass a suffix so plots/files are named uniquely
                self.filename_suffix = suffix
                self.run_process_chain(suffix=self.filename_suffix)
                
                # Clean up memory/cache for next window
                self.mseed_stream = None 
                self.dfdas = None
            else:
                print(f"Skipping window {i+1} due to loading error.")
    
        return False

    def load_event_data(self, event, mode='seismic', replace_DAS=False):
        """
        The main pipeline execution method. 
        """
        print(f'Working on event: {event}.')

        # 1. Initial setup
        self.dfdas = None
        self.dfdas_offsets = None
        self.mseed_stream = None
        self.event = event
        station = self.exec_params.get("station")
        
        # Accessing sources (ensure sources is available globally or passed in)
        self.acquisition_params["gps_status"] = self.sources.get(event, {}).get("GPS", False)
        
        # 2. Get initial trim times
        if mode != "time_lapse":
            # Standard logic using station files
            self.figsize, self.datetime_start, self.duration = self.calculate_trim_times(event,station)
        else:
            # Time Lapse Logic: 
            # self.datetime_start was already set by process_time_lapse loop.
            # We MUST ensure it is a native datetime object for simpledas compatibility.
            if hasattr(self.datetime_start, 'datetime'):
                self.datetime_start = self.datetime_start.datetime
            self.figsize = (10, 8) # Default for time-lapse

        # 3. Locate files and metadata
        self.locate_event_data(event)
        if not self.files:
            print(f"No files found for {event}. Execution aborted for this event.")
            return False 
        
        # 4. Handle Stale GPS Logic
        # We check if the GPS status changed after looking at the file headers
        if self.acquisition_params.get("gps_status") == True:
            if self.sources.get(event, {}).get("GPS") == False:
                print("GPS information was stale. Recalculating start times.")
                figsize, self.datetime_start, self.duration = self.calculate_trim_times(event, station)
                # Re-locate files with updated time
                self.locate_event_data(event)

        # 5. Check if we found anything
        if not self.files:
            print(f'No data available for event {event}. Skipping!')
            return False # Signal to the loop to continue

        # 6. Calculate Channels
        self.channels, min_ch, max_ch, ch_step = self.calculate_channels(event)

        # Exit early if channels array is empty
        if self.channels.size == 0:  # Use .size for NumPy arrays
            print(f"No valid channels found for event {event}. Exiting.")
            return False

        # 7. Final Load
        self.load_dfdas_file(event, replace_DAS=replace_DAS)

        if self.exec_params.get("verbose"):
            print(f"LOAD_EVENT_DATA: Event {event} loaded successfully with {len(self.channels)} channels.")
        
        return True # Success!

    def locate_event_data(self,event):

        events_path = self.path_params.get("events_path")
        experiment_path = self.path_params.get("experiment_path")
        # replace spaces with underscores in event name for filename
        hdf5_file = f'{events_path}/{event.replace(" ", "_")}.hdf5'

        self.event=event

        print(f'Finding files for event: {event}')
        print(f'Trying file list from {events_path}')

        if os.path.exists(hdf5_file):
            self.files=[hdf5_file]
            with h5py.File(self.files[0], "r") as hf_in:
                self.roi_step=hf_in['/demodSpec/roiDec'][0]
                self.first_chan = hf_in['/header/channels'][0]
                self.last_chan = hf_in['/header/channels'][-1]
                ap = self.acquisition_params 
                ap["roi_step"]=hf_in['/demodSpec/roiDec'][0]
                ap["gps_status"]=hf_in['/monitoring/Gps/gpsStatus'][()]
                ap["gauge_length"]=hf_in['/header/gaugeLength'][()]
                ap["dx"]=hf_in['/header/dx'][()]
                
                # Calculate overlap using the newly updated dictionary
                ap["overlap"] = ap["roi_step"] * ap["dx"] / ap["gauge_length"]
        else:
            print(f'Local HDF5 not found. Searching raw data in {experiment_path}')
            date_start= self.datetime_start.strftime("%Y%m%d")
            target_path=Path(experiment_path) / date_start
            if not target_path.exists():
                print(f"Skipping event {event} due to missing directories.")
                self.files = []
                return self.files, 0, 0, 0
           
           # If check passes, proceed to raw data metadata extraction 
            self.get_roi_step(experiment_path, self.datetime_start)

        return self.files, self.roi_step, self.first_chan, self.last_chan

    def get_roi_step(self, path, datetime_start):
        """
        Extracts ROI and acquisition metadata from DAS files.
        Updates self.acquisition_params directly.
        """
        # Call simpledas using parameters stored in the class
        files, chIndex, samples = simpledas.find_DAS_files(
            path,
            start=datetime_start,
            duration=datetime.timedelta(hours=1),
            show_header_info=self.exec_params.get("verbose"),
            load_file_from_start=True
        )

        if len(files) == 0:
            self.files = []
            self.roi_step = 0
            self.last_chan = 0
            self.first_chan = 0
        else:
            self.files = files
            with h5py.File(self.files[0], "r") as hf_in:
                # Update individual attributes
                self.first_chan = hf_in['/header/channels'][0]
                self.last_chan = hf_in['/header/channels'][-1]
                self.roi_step = hf_in['/demodSpec/roiDec'][0]

                # Update the shared acquisition dictionary
                ap = self.acquisition_params # Shorthand reference
                ap["roi_step"] = self.roi_step
                ap["gps_status"] = hf_in['/monitoring/Gps/gpsStatus'][()]
                ap["gauge_lenght"] = hf_in['/header/gaugeLength'][()]
                ap["dx"] = hf_in['/header/dx'][()]
                ap["original_sample_rate"] = 1.0 / hf_in['/header/dt'][()]
                ap["sensitivity"] = hf_in['/header/sensitivities'][0][0]
                ap["sensitivityUnits"] = hf_in['/header/sensitivityUnits'][0].decode('utf-8')
                
                # Calculate overlap
                ap["overlap"] = ap["roi_step"] * ap["dx"] / ap["gauge_lenght"]

        return self.files, self.roi_step, self.first_chan, self.last_chan

    def calculate_trim_times(self, event, station=None):
        """
        Calculates figure size, start time, and duration for the event.
        Stores results in self.figsize, self.datetime_start, and self.duration.
        """
        ap = self.acquisition_params
        ep = self.exec_params
        einfo = self.event_info
        if station is None:
            station = ep.get("station")


        if self.roi_step is None:
            logging.warning("ROI step is not set. Attempting to retrieve from loaded data.")
            self.roi_step= self.filemeta['demodSpec']['roiDec'][0]
            if self.roi_step is None:
                logging.warning("ROI step could not be retrieved from current dfdas data. Attempting to retrieve from files.")
                experiment_path = self.path_params.get("experiment_path")
                start_time = UTCDateTime(einfo.get("start_time")).datetime
                self.get_roi_step(experiment_path, start_time)
            if self.roi_step is None:
                logging.error("ROI step could not be determined. Please check the data files.")
                raise ValueError("ROI step is not set and could not be retrieved from files.")

        print(f'Processing event {event}, of the type {einfo.get("type")}')
        print('Calculating trim times and figure sizes...')

        # Determine clock offset based on GPS status stored in the class
        if ap.get("gps_status"):
            clock_offset = 0
        else:
            # Note: We assume clock_offset was passed into global_params 
            # and stored in self at __init__
            clock_offset = ep.get("clock_offset", 0) 

        # Logic for Seismic vs Other types
        if einfo.get("type", "").startswith('seismic'):
            self.figsize = (10, 3)
            # Assuming get_travel_times is an external helper
            arrivals_p, arrivals_s = get_travel_times(
                event, 
                self.get_current_params(), # Helper to package params if needed
                self.station_coordinates[station]["latitude"], 
                self.station_coordinates[station]["longitude"]
            )
            
            # Calculate start time with offset
            base_time = UTCDateTime(einfo.get("start_time"))
            trim_time = base_time + arrivals_p[0].time - einfo.get("window_in") + datetime.timedelta(seconds=clock_offset)
            print(f'Time offset correction of {clock_offset}s applied.')
        else:
            self.figsize = (10, 8)
            trim_time = UTCDateTime(einfo.get("start_time")) - einfo.get("window_in")
            print("Time offset not corrected.")

        # Store results in the class instance
        self.datetime_start = trim_time.datetime
        self.duration = datetime.timedelta(seconds=einfo.get("window_duration") + einfo.get("window_in"))

        print(f'Datetime start: {self.datetime_start}')
        return self.figsize, self.datetime_start, self.duration

    def trim_DAS_data(self, event):
        """
        Trims the DAS data to the specified time window.
        """

        window_in = pd.to_timedelta(self.event_info.get("window_in"), unit='s')
        window_duration = pd.to_timedelta(self.event_info.get("window_duration"), unit='s')
        start_time = pd.to_datetime(self.event_info.get("start_time"))

        self.dfdas = self.dfdas.truncate(before=start_time - window_in, after=start_time + window_duration)
        print('DAS data trimmed from ', start_time - window_in, ' to ', start_time + window_duration)
        print('with a start time of ', start_time)

        return self.dfdas

    def load_dfdas_file(self, event, replace_DAS=False):
        """
        Loads the actual DAS data into self.dfdas. 
        Uses the results from locate_event_data if they exist.
        """
        # Accessing class attributes
        ep = self.exec_params
        pp = self.path_params
        ap = self.acquisition_params
        
        datetime_start = self.datetime_start
        duration = self.duration
        channels = self.channels
        
        einfo = self.event_info
        event_type = einfo.get("type")
        event_mode = einfo.get("mode")

        is_special_type = event_type in ['spectrogram', 'whales', 'spectrogram_currents']
        #local_hdf5 = f'{pp.get("events_path")}/{event}.hdf5'
        # replace spaces with underscores in event name for filename
        local_hdf5 = f'{pp.get("events_path")}/{event.replace(" ", "_")}.hdf5'

        # 1. Logic to decide WHERE to load from
        if is_special_type or event_mode == "time_lapse" or replace_DAS:
            self.files, ch_idx, samples = simpledas.find_DAS_files(
                pp.get("experiment_path"), channels=channels, 
                start=datetime_start, duration=duration, load_file_from_start=False
            )
            self.dfdas = simpledas.load_DAS_files(self.files, ch_idx, samples, integrate=False)
        else:
            if os.path.exists(local_hdf5):
                try:
                    self.files = [local_hdf5]
                    # This avoids channel re-numbering
                    ch_idx, samples = simpledas.get_data_indexes(self.files[0])
                    print(self.files[:])
                    self.dfdas = simpledas.load_DAS_files(self.files, ch_idx, samples, integrate=False)
                except ValueError:
                    # Fallback to experiment path if local HDF5 fails
                    self.files, ch_idx, samples = simpledas.find_DAS_files(
                        pp.get("experiment_path"), channels=channels, 
                        start=datetime_start, duration=duration, load_file_from_start=False
                    )
                    self.dfdas = simpledas.load_DAS_files(self.files, ch_idx, samples, integrate=False)
            else:
                self.files, ch_idx, samples = simpledas.find_DAS_files(
                    pp.get("experiment_path"), channels=channels, 
                    start=datetime_start, duration=duration, load_file_from_start=False
                )
                self.dfdas = simpledas.load_DAS_files(self.files, ch_idx, samples, integrate=False)

        # 2. Sanity check using acquisition_params (updated by previous methods)
        if einfo.get("GPS") and ap.get('gps_status') != 0:
            print("Event defined as having locked GPS, but file header indicates otherwise.\nClock offset should be calculated.")

        if ep.get("verbose"):
            print(f'chIndex: {ch_idx}')
            print(f'columns: {self.dfdas.columns}')

        # 3. Store Metadata
        self.filemeta = simpledas.get_filemeta(self.files[0], 2)

        # 3.1 update self.roi_step in case it was not set before
        if self.roi_step is None:
            logging.warning("ROI step was not set upon loading. Retrieving from file metadata.")
            self.roi_step = self.filemeta['demodSpec']['roiDec'][0]
            ap["roi_step"] = self.roi_step

        # 3.2 Update self.first_chan and self.last_chan if they were not set before
        if self.first_chan is None:
            logging.warning("First channel was not set upon loading. Retrieving from file metadata.")
            self.first_chan = self.filemeta['header']['channels'][0]
            ap["first_chan"] = self.first_chan
        if self.last_chan is None:
            logging.warning("Last channel was not set upon loading. Retrieving from file metadata.")
            self.last_chan = self.filemeta['header']['channels'][-1]
            ap["last_chan"] = self.last_chan

        # 4. Save if required
        if ep.get("save_DFDAS") and event_mode != "time_lapse":
            os.makedirs(pp.get("events_path"), exist_ok=True)
            file_exists = os.path.exists(local_hdf5)
            
            if file_exists and not replace_DAS:
                print(f'Local HDF5 already exists for event {event}. Skipping save.')
            else:
                if file_exists:
                    print(f'Local HDF5 already exists for event {event}. Overwriting.')

                filename_out = simpledas.save_to_DAS_file(
                    self.dfdas, local_hdf5, "processed", **self.filemeta
                    )
                print(f'Window for event {event} saved to {filename_out[0]}')

        # 5. populate offsets to make it easier to track channels in the DataFrame. This is useful for plotting and analysis.
        self.dfdas_offsets = self.dfdas.columns.to_numpy()

        return self.files, self.dfdas, self.filemeta
    
#TODO: This method is not currently called in the chain, but it can be used after loading to check for gaps before proceeding with processing.
# It looks at the 't' column (time axis) in self.dfdas, calculates the time differences between consecutive samples, and identifies any gaps that exceed a certain threshold based on the expected sampling interval. If gaps are found, it returns a list of their indices and durations for further handling (e.g., interpolation, exclusion).
# The threshold for identifying gaps can be adjusted based on the specific requirements of the analysis and the characteristics of the data.
# Note: This method assumes that self.dfdas has a 't' column representing time, which is typically created during the decimation process. If the data is loaded without this column, the method will simply return an empty list and print a warning.
# 
    def check_for_gaps(self):
        """
        Checks for gaps in the time series data based on the 't' column.
        Returns a list of gap intervals if found.
        """
        if 't' not in self.dfdas.columns:
            print("No 't' column found in DataFrame for gap checking.")
            return []

        time_diffs = self.dfdas['t'].diff().dropna()
        expected_dt = self.acquisition_params.get("roi_step", 1) * self.acquisition_params.get("dx", 1) / self.acquisition_params.get("original_sample_rate", 1)
        
        # Define a threshold for identifying gaps (e.g., 1.5 times the expected dt)
        gap_threshold = expected_dt * 1.5
        gaps = time_diffs[time_diffs > gap_threshold]

        if not gaps.empty:
            print(f"Warning: Detected {len(gaps)} gaps in the data exceeding {gap_threshold} seconds.")
            return list(zip(gaps.index, gaps.values))
        else:
            print("No significant gaps detected in the data.")
            return []

    def save_dfdas_file(self, params=None):
        """
        Saves self.dfdas to an HDF5 file using self.filemeta.
        Automatically constructs the path if 'filename' is not in params.
        """
        if self.dfdas is None:
            print("No data in memory (self.dfdas) to save.")
            return

        # 1. Determine the filename
        # Check if 'filename' was passed in the process_chain dictionary
        p = params if isinstance(params, dict) else {}
        filename = p.get("filename")

        if not filename:
            # Fallback: Construct filename from class attributes
            pp = self.path_params
            event_safe = self.event.replace(" ", "_")
            # Uses the suffix we set in the time_lapse loop
            suffix = getattr(self, 'filename_suffix', "")
            
            filename = os.path.join(
                pp.get("events_path"), 
                f"{event_safe}{suffix}.hdf5"
            )

        # 2. Execute save using class attributes
        try:
            filename_out = simpledas.save_to_DAS_file(
                self.dfdas, 
                filename, 
                "processed", 
                **self.filemeta
            )
            print(f"Data successfully saved to: {filename_out[0]}")
        except Exception as e:
            print(f"Error saving DFDAS file: {e}")

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
            print("No valid offsets are available for plotting.")
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

    def extract_channel_subset(self, offset_range_key):
        """
        Extracts a subset of channels from self.dfdas based on provided indices.
        Updates self.dfdas in place.
        """

        logging.info(f"Extracting channel subset for offset range key: {offset_range_key}")

        #1. create the offsets array
        ranges = self.offset_ranges.get(offset_range_key, {})
        offsets = []
        for min_offset, max_offset, step_offset in zip(
            ranges.get("offset_min", []),
            ranges.get("offset_max", []),
            ranges.get("offset_step", [])
        ):
            offsets.extend(np.arange(min_offset, max_offset + 1, step_offset))

        #         
        if self.dfdas is None:
            print("No data loaded in self.dfdas. Cannot extract channels.")
            return

        if self.dfdas_offsets is None:
            self.dfdas_offsets=self.dfdas.columns.to_numpy()

        original_offsets=set(self.dfdas_offsets)

        subset_dfdas = self.dfdas.columns.intersection(offsets)
        if subset_dfdas.empty:
            print("No matching channels found in self.dfdas for the provided offsets.")
            return

        missing_channels= set(offsets) - set(subset_dfdas)
        if missing_channels:
            # create a map of missing channels to their corresponding offset ranges for better clarity
            plot_channel_map = self.plot_channel_map(offsets, original_offsets)

            print(f"Warning: {len(missing_channels)} requested channels were not found in self.dfdas, and were skipped.")
            print(f"Successfully extracted: {len(subset_dfdas)} channels.")
            logging.warning(f"The following requested channels were not found in self.dfdas: {sorted(list(missing_channels))}")
            logging.info(f"Successfully extracted: {len(subset_dfdas)} channels.")

        self.dfdas = self.dfdas[subset_dfdas]
        self.dfdas_offsets = self.dfdas.columns.to_numpy()

        return

    def calculate_channels(self, event, time_lapse=False):
        """
        Calculates the channels array based on ROI and file headers.
        Uses attributes stored during locate_event_data().
        """
        event_channel_ranges = self.event_info.get("channel_ranges")

        if not event_channel_ranges:
            print(f"WARNING: 'channel_ranges' not specified for event '{event}'. Skipping channel calculation.")
            self.channels = np.array([])
            return self.channels, [], [], [] #This will signal an error upstream

        # Get ranges from the global channel_ranges dict
        ranges = self.channel_ranges.get(event_channel_ranges, {})
        self.min_channel = ranges.get("channel_min", [])
        self.max_channel = ranges.get("channel_max", [])
        self.channel_step = ranges.get("channel_step", [])
        
        temp_channels = []
        
        # Use self.roi_step and self.last_chan stored from previous methods
        for min_ch, max_ch, step_ch in zip(self.min_channel, self.max_channel, self.channel_step):
            top = (max_ch + 1) * self.roi_step if (max_ch + 1) < self.last_chan else self.last_chan
            temp_channels.extend(np.arange(min_ch * self.roi_step, top, step_ch * self.roi_step))
        
        self.channels = np.array(temp_channels)
        
        # Adjustment for non-zero first channel
        if self.first_chan != 0 and self.first_chan < self.roi_step:
            adjustment = self.first_chan - self.roi_step
            print(f'Adjusting channels: {self.channels[0]} -> {self.channels[0] + adjustment}')
            self.channels = self.channels + adjustment
            
        return self.channels, self.min_channel, self.max_channel, self.channel_step

    def dfdas_to_csv(self, params=None):
        """
        Saves the decimated self.dfdas directly to CSV.
        Uses 'csv_out_path' and appends the timelapse suffix.
        """
        if self.dfdas is None:
            print("No data to save.")
            return

        pp = self.path_params

        # 1. Setup path
        csv_dir = pp.get("csv_out_path")
        if not os.path.exists(csv_dir):
            os.makedirs(csv_dir)

        # 2. Build filename
        event_safe = self.event.replace(" ", "_")
        filename = f"{event_safe}{self.filename_suffix}.csv"
        full_path = os.path.join(csv_dir, filename)

        # 3. Export
        print(f"»»» Saving CSV format to: {full_path}")
        # index=True includes the timestamps from the 't' axis created in decimation
        self.dfdas.to_csv(full_path, index=True, header=True)


    def dfdas_to_mseed(self,event,write=False,format=None,dType=None,decimate=False,filename=None, stationXML=False):
        """
        Converts the stored self.dfdas (Pandas) into an ObsPy Stream object.
        Uses class attributes for metadata and configuration.
        """

        if dType == 'int':
            print(f'Scaling data by a factor of 1x10^12 to match int32 range for MSEED output.\nNew data type is {dType}.')
        print(f'params: event:{event} write:{write} format:{format} dType:{dType} decimate:{decimate} filename:{filename}')
        if self.dfdas is None:
            print("No data loaded in self.dfdas. Run load_dfdas_file() first.")
            return None
        else:
            print("Converting to ObsPy Stream")

        # Short aliases for cleaner code
        ep, ap, pp, th = self.exec_params, self.acquisition_params, self.path_params, self.trace_header[self.network]

        self.stationxml_params["int32_resolution"]='picostrain'

        self.mseed_stream=Stream()
        for channel in self.dfdas.columns[:]:
            tr = Trace()
            tr.data = self.dfdas[channel].to_numpy().astype(float)  # Convert to dType. default is float.
            if dType=="int":
                tr.data = tr.data*1e6  # Scale to int32 range for MSEED. 
                self.stationxml_params["int32_resolution"]='picostrain'
                tr.data = tr.data.astype(np.int32) # Convert to int32 for MSEED
            tr.stats.starttime = UTCDateTime(self.dfdas.meta["time"])
            tr.stats.sampling_rate = 1/self.dfdas.meta["dt"]

            station_code = str(int(channel*ap.get("dx"))).zfill(5)
            if ep.get("encode_station_name"):
                tr.stats['decoded_station'] = station_code
                tr.stats.station = encode_station_name(th.get("station_prefix"),station_code,self.get_current_params())
            else:
                tr.stats.station = station_code
            tr.stats.network = th['network']
            tr.stats.location = th['location']
            tr.stats.mseed = th['mseed']
            tr.stats.channel = th['channel']
            tr.stats.distance = int(channel*self.dfdas.meta["dx"])
            lat, long, ele = get_coordinates(channel, self.get_current_params(),roi_channel=False)
            
            if ep.get("degrade_coordinates"):
                lat,long=degrade_coordinates(lat,long,precision=ep.get("precision"),error_metres=ep.get("error_metres"))
            
            tr.stats.coordinates = {'latitude': lat, 'longitude': long, 'elevation': ele}
            tr.stats.paz= AttribDict({'sensitivity': ap["sensitivity"]})
            # Append the trace to the Stream
            self.mseed_stream.append(tr)

        if decimate:
            self.mseed_stream=self.decimate_st(self.mseed_stream,self.get_current_params())

        if ep.get("stationXML") or stationXML:
            if ep.get("verbose"):
                print("Creating StationXML file for the current stream.")
                station_xml_creation(self.mseed_stream,self.get_current_params())   

        if write in ['stream', True]:
            mseed_filename = filename if filename else f'{event}_{self.mseed_stream[0].stats.sampling_rate}Hz.mseed'
            full_path = os.path.join(pp.get("MSEED_path"), mseed_filename)
            self.mseed_stream.write(full_path, format='MSEED')
            print(f'Full stream written to {full_path}')

        elif write == 'channels':
            for tr in self.mseed_stream:
                net = tr.stats.network  # Network code
                sta=tr.stats.station
                loc = tr.stats.location or ""  # Location code (or empty string if not provided)
                chan = tr.stats.channel  # Channel code
                qual= tr.stats.mseed # Data quality flag
                if ep.get("julian_date"):
                    formatted_time=f'{tr.stats.starttime.year:04d}.{tr.stats.starttime.julday:03d}'
                else:
                    formatted_time = tr.stats.starttime.strftime('%Y%m%dT%H%M%S')

                mseed_filename=f'{net}.{sta}.{loc}.{chan}.{qual}.{formatted_time}.mseed'
                if dType=="int":
                    mseed_path=os.path.join(pp.get("MSEED_path"),'int32',f'{tr.stats.starttime.year:04d}.{tr.stats.starttime.julday:03d}_int')
                else:
                    mseed_path=os.path.join(pp.get("MSEED_path"),f'{tr.stats.starttime.year:04d}.{tr.stats.starttime.julday:03d}')
            
                if not os.path.exists(mseed_path):
                    os.makedirs(mseed_path)
            
                full_path=os.path.join(mseed_path,mseed_filename)
                tr.write(full_path, format='MSEED',reclen=512)#,encoding='STEIM1')
                print(f'wrote file {full_path} julian_date={ep.get("julian_date")}')

        elif not write:
            print("Stream created in memory. No files written.")

        return self.mseed_stream
    
    def cross_correlate(self, reference=None, max_lag_seconds=10):
        """
        Cross-correlates all traces in self.dfdas and stores results in self.cross_correlation_results.
        """
        if self.dfdas is None:
            print("No data is loaded. Run load_dfdas_file() first.")
            return None
        else:
            print("Performing cross-correlation on self.dfdas.")
        
        if reference is None:
            num_channels = len(self.dfdas.columns)
            reference_idx= (num_channels - 1) // 2
        else:
            reference= int(reference/self.roi_step) # Convert from channel number to index based on roi_step
            if reference >= len(self.dfdas.columns) or reference < -len(self.dfdas.columns):  # the negative index allows counting from the end, but it should still be within bounds
                print(f"Reference channel index {reference} is out of bounds. Defaulting to middle channel.")
                num_channels = len(self.dfdas.columns)
                reference_idx = (num_channels - 1) // 2
            else:
                reference_idx = reference

        reference_channel = self.dfdas.columns[reference_idx]
        channel_distance = int(reference_channel * self.dfdas.meta['dx'])
        print(f"Using channel no. {reference_idx} as the reference channel for cross-correlation.")
        print(f"Distance for reference channel: {channel_distance} metres.")

        self.cross_correlation_results = {}

        reference_signal = self.dfdas[reference_channel].to_numpy()
        sampling_rate = 1.0 / self.dfdas.meta['dt']
        max_lag_samples = int(max_lag_seconds * sampling_rate) 
        lags = sps.correlation_lags(len(reference_signal), len(reference_signal), mode='full') / sampling_rate
        valid_indices = np.where(np.abs(lags) <= max_lag_seconds)
        trimmed_lags = lags[valid_indices]

        for channel in self.dfdas.columns:
            signal = self.dfdas[channel].to_numpy()
            corr = sps.correlate(signal, reference_signal, mode='full')
            # Limit to max_lag_seconds
            self.cross_correlation_results[channel] = {
                'correlation': corr[valid_indices],
                'lags': trimmed_lags
            }
        print(f"Cross-correlation completed with reference channel {reference_channel}.")
        return self.cross_correlation_results

    def plot_cross_correlation_heatmap(self, vmin=None, vmax=None, cmap='RdBu_r'):
        """
        Plots all cross-correlation results as a 2D heatmap/waterfall plot.
        """
        if not hasattr(self, 'cross_correlation_results') or not self.cross_correlation_results:
            print("No cross-correlation results found. Run cross_correlate() first.")
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
            vmax = np.percentile(np.abs(corr_matrix), 98)  # Clip top 2% outliers for better contrast
            vmin = np.percentile(np.abs(corr_matrix), 2)   # Clip bottom 2% outliers for better contrast

        im = plt.imshow(corr_matrix, aspect='auto', extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)

        plt.colorbar(im, label='Correlation Coefficient')
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
            print("No cross-correlation results found.")
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
    
    
    def decimate_st(self):
        ep, fp, sp = self.exec_params , self.freq_params, self.stationxml_params

        st=self.mseed_stream

        for decimation in fp.get("decimation_factor"):
            print(f'Decimation factor: {decimation}')
            new_nyquist = st[0].stats.sampling_rate / 2.0 / decimation
            b,a=get_cheby2_coeffs(st[0].stats.sampling_rate * 0.5,new_nyquist)
            sp["input_sampling_rate"].append(st[0].stats.sampling_rate)
            sp["filter_a_coefficients"].append(a)
            sp["filter_b_coefficients"].append(b)
            for tr in st:
                zerophase_chebyshev_lowpass_filter(tr, new_nyquist)
            st.decimate(factor=decimation, no_filter=True)
            if ep.get("verbose"):
                print(st[0].stats.sampling_rate)
        return st

    def reduce_velocity_st(self, global_params, vred=5.0, offset=None):
        ep, fp = global_params['exec_params'] , global_params['freq_params']

        st=self.mseed_stream.copy()
        if offset is not None:
            print(f'An offset of {offset} meters with be added to every distance.')

        for tr in st:
            # apply a distance shift to match a specific feature
            tr.stats.distance += offset

            # Calculate the time shift based on distance
            # distance is usually in meters in ObsPy, so convert vred to m/s
            shift = tr.stats.distance / (vred * 1000)

            # Update the starttime to shift the trace relative to origin
            tr.stats.starttime -= shift

        self.mseed_reduced_velocity = st  # Store the reduced velocity stream in the class for later use
        return st

    def decimate_dfdas(self, factors=None):
        """
        Decimates self.dfdas based on the decimation factors provided in configuration.
        Updates self.dfdas and its metadata in-place.
        """

        fp, sp = self.freq_params, self.stationxml_params

        if self.dfdas is None:
            print("Data not found in memory. Attempting to load before decimation...")
            return False
        
        # If factors are passed via process_chain, use them; 
        # otherwise, fall back to global_params
        if factors is None:
            factors = fp.get("decimation_factor", [])
        else:
            fp["decimation_factor"] = factors  # Update global_params if factors are passed directly
            print(f"Default decimation factors updated to: {factors}")
        
        for decimation in factors:
            print(f'Decimating at a factor of {decimation}')
            
            # update stationXML params for each decimation step
            original_sampling_frequency = 1/self.dfdas.meta['dt']
            new_nyquist = original_sampling_frequency * 0.5 / decimation
            b,a=get_cheby2_coeffs(original_sampling_frequency,new_nyquist)
            sp["input_sampling_rate"].append(original_sampling_frequency)
            sp["filter_a_coefficients"].append(a)
            sp["filter_b_coefficients"].append(b)
            
            # 1. Capture current dimensions
            nt_in = self.dfdas.shape[0]
            
            # 2. Perform Resampling
            # resample_poly includes an anti-aliasing filter automatically
            sig_decimated = sps.resample_poly(
                self.dfdas, up=1, down=decimation, axis=0, padtype='edge'
            )
            
            # 3. Calculate new timing metadata
            nt_out = sig_decimated.shape[0]
            dt_out = self.dfdas.meta['dt'] * (nt_in / nt_out)
            tstart = self.dfdas.meta['time'] # + datetime.timedelta(seconds=0)
            
            # Use your helper to create the new index
            t = simpledas.create_time_axis(tstart, nt_out, dt_out)
            
            # 4. Update metadata and re-wrap as a DASDataFrame
            meta_out = self.dfdas.meta.copy()
            meta_out.update(dt=dt_out, time=tstart)
            
            # Overwrite the class attribute with the new decimated version
            self.dfdas = simpledas.DASDataFrame(
                sig_decimated, 
                index=t, 
                columns=self.dfdas.columns, 
                meta=meta_out
            )
            
            print(f"New sampling rate: {1/dt_out} Hz")
            # CRITICAL: Clear the cached ObsPy stream because the underlying 
            # data (sampling rate/length) has changed.
            self.mseed_stream = None 

        return self.dfdas

    def plot_spectrogram(self, freq_type="no_filter", suffix=""):
        ep = self.exec_params
        pp = self.path_params
        sfp = self.exec_params.get('spectrogram_freq_params', {})

        if self.mseed_stream is None:
            print("No ObsPy Stream found in memory. Attempting to create from self.dfdas...")
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
        print("saving figure to ",filename)
        distance_km = trace.stats.distance / 1000
        fig.suptitle(f'Spectrogram for channel at {distance_km:.3f} km', y=1.0)
        fig.savefig(filename+'.pdf',dpi=600, format='pdf', bbox_inches='tight', pad_inches=0.1)
        fig.savefig(filename+'.png',dpi=300, format='png')
        plt.show()
        plt.close()  

        return

    def calculate_hilbert_envelope(self, params=None):
        """
        Calculates and saves the signal envelope plot.
        Supports overrides for cmap, vmin, and vmax via params.
        """
        if self.dfdas is None:
            print("No data in self.dfdas to plot.")
            return

        # 1. Setup local references and shorthand
        ep = self.exec_params
        pp = self.path_params
        ap = self.acquisition_params
        
        # Ensure params is a dict even if None/list was passed
        p = params if isinstance(params, dict) else {}

        # 2. Determine Plotting Parameters (Override logic)
        # It checks 'p' first, then 'ep', then provides a hard default
        cmap = p.get("envelope_cmap", ep.get("envelope_cmap", "viridis"))
        vmin = p.get("envelope_vmin", ep.get("envelope_vmin", .7))
        vmax = p.get("envelope_vmax", ep.get("envelope_vmax", 2))

        print(f"»»» Running envelope (Cmap: {cmap}, Vmin: {vmin}, Vmax: {vmax})")
        # 1. Attempt to extract the real part; if it's a Stream object, this will fail.
        try:
            data = self.dfdas.real
        except AttributeError:
        # 2. If it fails, we fall back to using 'data' directly
            data = self.dfdas

        en = abs(sps.hilbert(data,axis=1))
        mdn = np.median(en, axis=1, keepdims=True)
        envelope = en / abs(mdn)

        return envelope

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
        print(f"Envelope plot saved as {filename}.pdf")

    def check_data_access(self, experiment_path, datetime_start):
        """
        Verifies that the experiment path and the date-specific 
        subfolder exist before attempting to read files.
        """
        print(f"Checking data access in {experiment_path} for date {datetime_start.strftime('%Y-%m-%d')}")
        # 1. Check base path
        if not os.path.exists(experiment_path):
            print(f"ERROR: Base path does not exist: {experiment_path}")
            return False

        # 2. Check date-specific subfolder (SimpleDAS standard: path/YYYYMMDD)
        date_str = datetime_start.strftime('%Y%m%d')
        date_path = os.path.join(experiment_path, date_str)
        
        if not os.path.exists(date_path):
            print(f"ERROR: Date folder missing: {date_path}")
            print("Check if the drive is mounted or if the event date is correct.")
            return False

        # 3. Check for data type subfolder (e.g., 'dphi')
        dphi_path = os.path.join(date_path, 'dphi')
        if not os.path.exists(dphi_path):
             print(f"WARNING: 'dphi' folder not found in {date_path}. SimpleDAS might fail.")
             return False
        
        return True



    
    # ============================================================
    # PLOT METHODS
    # ============================================================

    def plot_seismic_section(self, channel_step=1, scale_factor=1.0, start_offset_km=30.8,
                              end_offset_km=None, dx_meters=12.5, max_time_s=45.0,
                              station=None, event=None, suffix="", use_agc=True,
                              agc_window_s=2.0, agc_clip=4.0, **kwargs):
        """
        Record section with wiggle traces, band-passed 0.5-20 Hz, with P/S/T phases.
        """
        if self.dfdas is None:
            print("No data in memory to plot.")
            return

        event = event or self.event
        print(f"»»» Generating Record Section for {event}...")

        data = self.dfdas.to_numpy()
        dt = self.dfdas.meta['dt']
        fs = 1.0 / dt
        window_start_utc = UTCDateTime(self.dfdas.meta['time'])
        time_axis = (self.dfdas.index - self.dfdas.index[0]).total_seconds()

        dx_native = self.dfdas.meta.get('dx', 1.0213001907746815)
        absolute_offsets_km = (self.dfdas.columns.to_numpy() * dx_native) / 1000.0

        if end_offset_km is None:
            end_offset_km = start_offset_km + 21.5

        print(f"-> Diagnostics: dx={dx_native:.4f} m/channel | dt={dt:.5f}s (fs={fs:.1f}Hz) | "
              f"available channel coverage: {absolute_offsets_km.min():.2f} - {absolute_offsets_km.max():.2f} km | "
              f"data window: {time_axis.min():.1f} - {time_axis.max():.1f} s | "
              f"window start (UTC): {window_start_utc}")

        spatial_mask = (absolute_offsets_km >= start_offset_km) & (absolute_offsets_km <= end_offset_km)
        sub_data = data[:, spatial_mask]
        offsets_km = absolute_offsets_km[spatial_mask]

        if sub_data.shape[1] == 0:
            print("ERROR: no channels in requested offset range.")
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
        print(f"-> Plotting {final_data.shape[1]} traces (target spacing: {dx_meters} m).")

        fig, ax = plt.subplots(figsize=(15, 8.5))
        max_deflection_km = trace_spacing_km * 0.45 * scale_factor
        agc_win_samples = max(1, int(agc_window_s * fs))
        agc_kernel = np.ones(agc_win_samples) / agc_win_samples

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

            ax.plot(final_offsets[i] + traccia_km, plot_time, color='#2c3e50', linewidth=0.4, alpha=0.8)

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
                        print(f"WARNING: phase {label} at t={t_val} is outside the plotted window "
                              f"({plot_time.min():.1f}-{plot_time.max():.1f}s).")
            except Exception as e:
                print(f"Unable to calculate theoretical phases: {e}")

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
        print(f"Figure saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_seismic_2d(self, channel_step=1, start_offset_km=30.8,
                        end_offset_km=50.0, dx_meters=12.5, max_time_s=45.0,
                        station=None, event=None, suffix="", **kwargs):
        """
        2D version of seismic plot using imshow.
        """

        if self.dfdas is None:
            print("No data in memory to plot.")
            return

        event = event or self.event
        print(f"»»» Generating 2D Record Section for {event}...")

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
            print(f"ERROR: no channels in range {start_offset_km}-{end_offset_km} km.")
            print(f"Available range: {absolute_offsets_km.min():.2f} - {absolute_offsets_km.max():.2f} km")
            return

        time_mask = time_axis <= max_time_s
        sub_data = sub_data[time_mask, :]
        plot_time = time_axis[time_mask]
        print(f"-> Data: {sub_data.shape[1]} channels, {sub_data.shape[0]} time samples")

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
                print(f"Unable to calculate theoretical phases: {e}")

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
        print(f"2D plot saved to: {filename}")
        plt.show()
        plt.close(fig)

    def plot_seismic_zoom(self, start_offset_km=30.8, end_offset_km=50.0,
                          t_min=5, t_max=30, event=None, suffix="", **kwargs):
        """
        Zoomed version of the seismic plot to view P and S arrivals.
        """

        if self.dfdas is None:
            print("No data in memory to plot.")
            return

        event = event or self.event
        print(f"»»» Generating ZOOM Record Section for {event}...")

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
                print(f"Error computing phases: {e}")

        plt.tight_layout()
        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        filename = os.path.join(save_path, f"{event}_record_section_zoom.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Zoom plot saved to: {filename}")
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
            print("No data in memory to plot.")
            return

        event = event or self.event
        print(f"»»» Generating ENHANCED Record Section for {event}...")

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
            print("ERROR: no channels in the requested range.")
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

        for i in range(sub_data.shape[1]):
            trace = sub_data[:, i]
            trace = sps.detrend(trace)
            trace = trace - np.median(trace)
            ax1.plot(offsets_km[i] + trace * max_deflection / 3.0, 
                    plot_time, 'k-', linewidth=0.3, alpha=0.7)

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
                print(f"Error computing phases: {e}")

        plt.suptitle(f"{event} - Seismic Record Section (M2.7 Desertas)", fontsize=14, fontweight='bold')
        plt.tight_layout()

        pp = self.path_params
        save_path = os.path.join(pp.get("figures_path", "./"), "Seismic_Plots")
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        filename = os.path.join(save_path, f"{event}_enhanced_section.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Enhanced plot saved to: {filename}")
        plt.show()
        plt.close(fig)

    def download_waveform_data(self, station_code=None, starttime=None, endtime=None,
                               network=None, location='*', channel='HHZ', eida_url=None,
                               save_stream=True):
        """
        Download waveform data for a land station from a hardcoded EIDA FDSN node.

        Prompts the user for missing parameters (station code, start/end times) when
        they are not provided. Times should be ISO-like strings parseable by ObsPy
        (e.g. '2020-01-01T00:00:00'). The downloaded Stream is returned and also
        stored on the instance under `self.station_streams[station_code]` if
        `save_stream` is True.

        Hardcoded EIDA base URL may be changed later; currently set to a default
        value compatible with EIDA webservices.
        """
        try:
            from obspy.clients.fdsn import Client
        except Exception:
            print("ObsPy FDSN client not available. Install ObsPy to enable downloads.")
            return None

        # Default EIDA node (can be changed later)
        if eida_url is None:
            eida_url = 'https://webdc.eu'  # hardcoded EIDA-compatible node

        # Ask for station code if not provided
        if station_code is None:
            station_code = input("Enter station code (e.g. ABCD): ").strip()
            if not station_code:
                print("No station code provided. Aborting download.")
                return None

        # Ask for network if not provided
        if network is None:
            net_in = input("Enter network code (or press Enter to use '*'): ").strip()
            network = net_in if net_in else "*"

        # Ask for start and end times if not provided
        if starttime is None:
            s = input("Start time (YYYY-MM-DDTHH:MM:SS): ").strip()
            try:
                starttime = UTCDateTime(s)
            except Exception as e:
                print(f"Invalid start time: {e}")
                return None
        else:
            starttime = UTCDateTime(starttime)

        if endtime is None:
            e = input("End time (YYYY-MM-DDTHH:MM:SS): ").strip()
            try:
                endtime = UTCDateTime(e)
            except Exception as ex:
                print(f"Invalid end time: {ex}")
                return None
        else:
            endtime = UTCDateTime(endtime)

        # Instantiate client and fetch waveforms
        try:
            client = Client(base_url=eida_url)
            print(f"Connecting to EIDA node: {eida_url}")
            st = client.get_waveforms(network, station_code, location, channel, starttime, endtime)
            print(f"Downloaded {len(st)} trace(s) for {station_code} from {starttime} to {endtime}.")
        except Exception as err:
            print(f"Error downloading data from EIDA: {err}")
            return None

        # Store stream on instance
        if save_stream:
            if not hasattr(self, 'station_streams') or self.station_streams is None:
                self.station_streams = {}
            self.station_streams[station_code] = st

        return st

    def fk_filter(self, fs, vmin=1400,vmax=99999):
        """
        Apply FK filter to the data
        """

        print(f"Applying FK filter with vmin={vmin}, vmax={vmax}.")

        # Get the Discrete Fourier Transform sample frequencies...
        dt = self.dfdas.meta['dt']
        dx = self.dfdas.meta['dx']
        nnx = self.dfdas.meta["dimensionSizes"][1]
        nns = self.dfdas.meta["dimensionSizes"][0]
        f = np.fft.fftfreq(nns,d=dt)
        k = np.fft.fftfreq(nnx,d=dx)

        f = np.fft.fftshift(f)
        k = np.fft.fftshift(k)

#        if self.mseed_stream is None:
#            self.dfdas_to_mseed(self.event, write=False, format=None, dType=float, decimate=False)
        
        # Compute the 2-dimensional discrete Fourier Transform...
        # Compute the 2D FFT and center the zero-frequency/wavenumber
        ft2s =np.fft.fftshift(
                              np.fft.fft2(self.dfdas) # check if it should be transposed or not
                              )

        # Return a list of coordinate matrices from coordinate vectors...
        F,K = np.meshgrid(f,k)
        del f, k

        # Change null values to very small ones to be able to divide afterwards.
        K = np.where(K==0, 1e-10, K)
        C = F/K

        # Create filter mask...
        # Create a mask of zeros in the same shape as C.
        filt = np.zeros(C.shape)
        # Create a mask to "mute" entries by replacing zeros with ones where C is in the desired range.
        condition = np.logical_and(np.abs(C) >= vmin, np.abs(C) <= vmax)
        filt[condition] = 1.

        # Blur mask.
        filt_blur = ndimage.gaussian_filter(filt, sigma=10)

        # Apply mask...
        ft2_filt = ft2s * filt_blur
        # del fts2
        ft2s_filt = np.fft.ifftshift(ft2_filt)
        del ft2_filt  # to save memory
        # Invert FFT.
        self.dfdas = np.fft.ifft2(ft2s_filt)

        return self.mseed_stream


    def plot_record_sections(self, dataset=None, event=None, suffix=""):
        # Use class attributes for configuration to integrate with DASProcessor
        ep, pp, fp = self.exec_params, self.path_params, self.freq_params
        station_coordinates = self.station_coordinates

        event_info = self.sources.get(event, {})
    
        print("Plotting record section.")
        result=subprocess.run(["bash","./plot_GMTmap.sh",
                                      str(event_info.get("longitude")),
                                      str(event_info.get("latitude")),
                                      f'"{str(event_info.get("magnitude"))} {event_info.get("depth")}km"',
                                      f'{str(event)}_map.ps'],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
        filename=f'{pp.get("figures_path")}/OptoDAS_{event}_record-section'
    
        if ep.get("verbose"):
            # Print the standard output and standard error
            print("Standard Output:")
            print(result.stdout)
            print("\nStandard Error:")
            print(result.stderr)
            
            # Check the return code
            if result.returncode == 0:
                print("Bash script executed successfully.")
            else:
                print(f"Bash script returned an error (return code {result.returncode}).")
    
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
        if event_info.get("local"):
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
        time_delta=event_info.get("window_in")+event_info.get("end_time")-event_info.get("start_time")
        
        pdiff_time=UTCDateTime(event_info.get("time"))+arrivals_p[0].time
        start_trim_time=pdiff_time-datetime.timedelta(seconds=event_info.get("window_in"))
        end_trim_time=start_trim_time+datetime.timedelta(seconds=time_delta)
    
        date_format = DateFormatter("%d-%m-%Y %H:%M:%S")
    
    
    
        # Reference offset relative to event origin time
        p_time = arrivals_p[0].time - event_info.get('window_in')
        event_time = UTCDateTime(event_info.get('time'))
    
        for arr in arrivals_p:
            # Keep time as float offset (seconds)
            arr.time = round(arr.time, 3)
            # Print absolute arrival time
            print(arr.name, event_time + arr.time)
            # Normalise relative arrival time
            arr.time -= p_time
            # Print updated absolute arrival time
            print(arr.name, event_time + arr.time)
    
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
                st = self.dfdas_to_mseed(event, write=False, decimate=ep.get("decimate_obspy"))
        elif isinstance(dataset, Stream):
            st = dataset
        else:
            # dataset is likely a pandas DataFrame (dfdas). Temporarily assign and convert.
            orig_dfdas = getattr(self, 'dfdas', None)
            self.dfdas = dataset
            restore_dfdas = True
            st = self.dfdas_to_mseed(event, write=False, decimate=ep.get("decimate_obspy"))

        if restore_dfdas:
            self.dfdas = orig_dfdas

        # If conversion didn't produce a Stream, bail out
        if st is None:
            print("No Stream available for plotting.")
            return

        st.filter("bandpass", freqmin=min_freq, freqmax=max_freq, corners=4, zerophase=True)
        
        # decimate DAS for plotting
        # Loop through each trace in the Stream
        # Create a new Stream object containing only every other trace
        every_other_trace = []
    
        if event_info.get("local"):
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
            print(f'Getting channels {station_coordinates.get(station,{}).get("channels")} from {ep.get("land_station")}')
            if "MSEED_file" in self.sources.get(event,{}):
                st_land=read(os.path.join(pp.get("MSEED_path", ""),self.sources.get(event, {}).get("MSEED_file", "")))
            else:
                st_land=read(pp.get("mseed_land_stations"))
            st_event=st_land.select(channel=station_coordinates.get(station,{}).get("channels"),station=ep.get("land_station"))
            st_event.trim(starttime=start_trim_time, endtime=end_trim_time)
        
        if len(st_event)>0:
            if start_trim_time<st_event[0].stats.starttime-0.01:
                print(f'{st_event[0].stats.station}.{st_event[0].stats.channel} has not enough data before the event.')
                print(f'{st_event[0].stats.starttime}>{start_trim_time} diff={start_trim_time-st_event[0].stats.starttime}')
                # window_correction=start_trim_time-st_event[0].stats.starttime
                # sources.get(event, {})["window_in"] = sources.get(event, {}).get("window_in")+window_correction
                # re_trim=True
            if end_trim_time<st_event[0].stats.endtime-0.01:
                print(f'{st_event[0].stats.station}.{st_event[0].stats.channel} has not enough data after the event.')
                print(f'{st_event[0].stats.endtime}>{end_trim_time} diff={end_trim_time-st_event[0].stats.endtime}')
                # window_correction=end_trim_time-st_event[0].stats.endtime
                # sources.get(event, {})["window_out"] = sources.get(event, {}).get("window_out")+window_correction
            # re_trim=True
        else:
            print(f'The selected channel ({station_coordinates.get(station,{}).get("channels")}) doesn\'t exist in {os.path.join(pp.get("MSEED_path", ""),self.sources.get(event, {}).get("MSEED_file", ""))}.')
            # Try to download land station data from a configured EIDA node if available
            try:
                print("Attempting to download land station data from configured EIDA node...")
                st_dl = self.download_waveform_data(station_code=station, starttime=start_trim_time, endtime=end_trim_time)
                if st_dl is not None and len(st_dl) > 0:
                    # Select requested channels and trim to window
                    try:
                        st_event = st_dl.select(channel=station_coordinates.get(station, {}).get("channels"), station=station)
                        st_event.trim(starttime=start_trim_time, endtime=end_trim_time)
                        print(f"Successfully downloaded and trimmed data for station {station}.")
                    except Exception as e:
                        print(f"Downloaded data present but required channels not found: {e}")
                else:
                    print("No data returned from EIDA for the requested interval/station.")
            except Exception as e:
                print(f"EIDA download attempt failed: {e}")
    
        if len(st_event)>0:
            unique_phases=[]
            st_event.detrend()
            st_event.taper(max_percentage=0.05, type='hann')
            st_event.filter("bandpass", freqmin=fp.get("land_station_min_freq"), freqmax=fp.get("land_station_max_freq"), corners=8, zerophase=True)
            print(f'Land stations bandpass filtered between {fp.get("land_station_min_freq")} and {fp.get("land_station_max_freq")} Hz')
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
                print(f'Adding {event_channel.stats.station}-{event_channel.stats.channel}')
                if event_channel.stats.station=="DAS":
                    ax1[idx].set_xlim(-ylim_stack, ylim_stack)
                    ax1[idx].set_title(f'{event_channel.stats.station}-{event_channel.stats.location}')
                else:
                    ax1[idx].set_xlim(-ylim, ylim)
                    ax1[idx].set_title(f'{event_channel.stats.station}-{event_channel.stats.channel}')
                ax1[idx].set_ylim(time_delta,0)
            
            label_space = 1 if self.sources.get(event, {}).get("local") else 10
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
                        if ep.get("verbose"):
                            print("UNIQUE ARRIVAL",arr.name,arr.time,plot_time,UTCDateTime(self.sources.get(event, {}).get("time"))+arr.time)
                        ax2[0].axhline(y=plot_time, color='red')
                        ax2[0].annotate(arr.name, xy=(ylim*0.9, plot_time), xytext=(-100, 0),
                                    textcoords='offset pixels')
            if ep.get("verbose"):
                print("Unique arrivals",unique_phases)
            
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
    
            pdiff_time=UTCDateTime(event_info.get("start_time"))+arrivals[0].time+datetime.timedelta(seconds=ep.get("clock_offset"))
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
        label_space = 1 if self.sources.get(event, {}).get("local") else 10
        if ep.get("verbose"): 
            print("LABEL SPACE:",self.sources.get(event, {}).get("local"),label_space)
    
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
    
        print(f'Interspace spacing={st_every_other[1].stats.distance-st_every_other[0].stats.distance} from {st_every_other[0].stats.distance} to {st_every_other[-1].stats.distance}. {len(st_every_other)} traces.')
        print("Finished plotting record section.")

def get_travel_times(event,global_params,rec_lat,rec_lon):
    """
    Compute theoretical P and S arrival times for an event at a receiver location.

    Parameters:
    - event: key identifying the source event in the `sources` dictionary.
    - global_params: dict containing at least 'exec_params' and 'sources'.
    - rec_lat, rec_lon: receiver latitude and longitude in degrees.

    Returns:
    - arrivals_p, arrivals_s: lists of TauP arrival objects for P and S phases.
    """
    ep = global_params['exec_params']
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
        
    if ep.get("verbose"):
        print(arrivals_p,arrivals_s)

    return arrivals_p, arrivals_s

def get_coordinates(channel, global_params,roi_channel=False):
    """
    Get coordinates (latitude, longitude, elevation) based on the provided channel

    Parameters:
    - channel (int): The channel number for which coordinates are needed.
    - coords_path (str): The path to the CSV file.
    - coords_DAS (str): The name of the CSV file.

    Returns:
    Tuple (float, float, float): A tuple containing latitude, longitude, and elevation.
    """
    ep, pp = global_params['exec_params'], global_params['path_params']
    
    if ep.get("verbose", False) == 2:
        print(f"Getting coordinates for channel {channel} from {pp['coords_path']}{pp['coords_DAS']}")
        
    df = pd.read_csv(os.path.join(pp['coords_path'],pp['coords_DAS']))

    if roi_channel:
        filtered_df = df.loc[df['ROI_CHANNEL'] == channel]
    else:
        filtered_df = df.loc[df['CHANNEL'] == channel]

   # TODO: replace hardcoded coordinates with a default value or raise an error if no match is found 
    if len(filtered_df) == 0:
        lat=32.64441
        lon=-16.94396
        ele=20
    else:
        lat=filtered_df['LATITUDE'].values[0]
        lon=filtered_df['LONGITUDE'].values[0]
        ele=round(filtered_df['ALTITUDE'].values[0],0)
    
    return lat,lon,ele

def encode_station_name(prefix,distance,global_params):
    """
    Encode a station name using a prefix and a numeric distance.

    Parameters:
    - prefix: string prefix for the station code.
    - distance: numeric distance (metres) used to generate the code suffix.
    - global_params: dictionary containing 'exec_params' for verbosity control.

    Returns:
    - code: formatted station code string.
    """
    ep=global_params.get("exec_params")   

    distance_m=int(distance)
    code_length = 5 - len(prefix)  # Calculate dynamic length
    factor=10 ** (5-code_length)
    converted_distance = distance_m // factor
    code = f"{prefix}{converted_distance:0{code_length}d}"

    if ep.get("verbose"):
        print(f'Encoded station: name: {code} - distance={distance} code_lenght={code_length} division_factor={factor}')
    return code


def station_xml_creation(st,global_params):
    """
    Create a StationXML file from an ObsPy Stream `st` using configuration in `global_params`.

    Parameters:
    - st: ObsPy Stream containing station traces with coordinates and metadata.
    - global_params: dictionary with acquisition, path, exec, freq, trace header and stationxml params.

    Side effects:
    - Writes StationXML to `pp.get("MSEED_path")` and may raise KeyError if coordinates missing.
    """
    ap, pp, ep  = global_params['acquisition_params'] , global_params['path_params'], global_params['exec_params']
    fp, th, sp = global_params['freq_params'], global_params['trace_header'], global_params['stationxml_params']
    
    #get filter params
    a,b=get_cheby2_coeffs(st[0].stats.sampling_rate*0.5,50)
    for tr in st:
        for coord in ["latitude","longitude","elevation"]:
            if coord not in tr.stats.coordinates:
                print(f'Station {tr.stats.station} has no {coord} field!\nStationXML file not created.')
                raise KeyError(f"Missing '{coord}' in tr.stats.coordinates")
                return

    response_stages=[]
    dummy_response_stage=[]
    
    response_stages=[PolesZerosResponseStage(stage_sequence_number=1,
                                             name="Gain",
                                             stage_gain=1,
                                             stage_gain_frequency=1.0,
                                             input_units="counts",
                                             output_units="counts",
                                             input_units_description=ap["sensitivityUnits"],
                                             output_units_description=ap["sensitivityUnits"],
                                             normalization_factor=1,
                                             normalization_frequency=0,
                                             pz_transfer_function_type="LAPLACE (RADIANS/SECOND)",
                                             zeros=[],
                                             poles=[])
                    ]

    for i in range(len(fp.get("decimation_factor"))):
    # Extract coefficients and parameters for the current stage
        decimation_factor = fp["decimation_factor"][i]
        input_sampling_rate = sp["input_sampling_rate"][i]
        a_coefficients = sp["filter_a_coefficients"][i].tolist()
        b_coefficients = sp["filter_b_coefficients"][i].tolist()

        # Create a CoefficientsTypeResponseStage object
        dummy_response_stage = CoefficientsTypeResponseStage(
            stage_sequence_number=i + 2,  # Stage sequence starts at 1
            name=f"Low-pass Chebishev filter stage + Decimation step {i + 1}. Decimation factor: {decimation_factor}",
            cf_transfer_function_type="DIGITAL",
            stage_gain_frequency=1.0,  # Frequency at which the gain is defined
            stage_gain=1,
            input_units="counts",
            output_units="counts",
            input_units_description=ap["sensitivityUnits"],
            output_units_description=ap["sensitivityUnits"],
            decimation_factor=decimation_factor,  # Set decimation factor
            decimation_input_sample_rate=input_sampling_rate,  # Input sample rate
            decimation_offset=0,
            decimation_delay=0,
            decimation_correction=0,
            denominator=a_coefficients,  # Denominator (a) coefficients
            numerator=b_coefficients,  # Numerator (b) coefficients
        )
        response_stages.append(dummy_response_stage)
    
    last_stage=len(response_stages)+1
    dummy_response_stage=CoefficientsTypeResponseStage(stage_sequence_number=last_stage,
                                  name="Sensitivity",
                                  cf_transfer_function_type="DIGITAL",
                                  numerator=[],
                                  denominator=[],
                                  decimation_factor=1,
                                  decimation_correction=0,
                                  decimation_delay=0,
                                  decimation_offset=0,
                                  decimation_input_sample_rate=ap.get("original_sample_rate")/math.prod(fp["decimation_factor"]),
                                  stage_gain=ap["sensitivity"],
                                  stage_gain_frequency=1.0,
                                  input_units="counts",
                                  output_units="counts",
                                  input_units_description=ap["sensitivityUnits"],
                                  output_units_description=ap["sensitivityUnits"],
                                  )

    response_stages.append(dummy_response_stage)

    print(sp.get("int32_resolution"))

    if sp.get("int32_resolution")== 'picostrain':
        picostrain_stage = CoefficientsTypeResponseStage(
                                 stage_sequence_number=last_stage + 1,
                                 name="Picostrain to Counts Resolution Mapping",
                                 cf_transfer_function_type="DIGITAL",
                                 stage_gain=1e12, 
                                 stage_gain_frequency=1.0,
                                 input_units="counts",
                                 output_units="PICOSTRAIN",
                                 input_units_description="Digital counts mapped to int32",
                                 output_units_description="Picostrain units (1e-12 strain)",
                                 numerator=[],
                                 denominator=[],
                                 decimation_factor=1,
                                 decimation_offset=0,
                                 decimation_delay=0,
                                 decimation_correction=0,
                                 decimation_input_sample_rate=ap.get("original_sample_rate")/math.prod(fp["decimation_factor"])
                )
        response_stages.append(picostrain_stage)
    
    if sp.get("int32_resolution")== 'microstrain':
        microstrain_stage = CoefficientsTypeResponseStage(
                                 stage_sequence_number=last_stage + 1,
                                 name="Microstrain to Counts Resolution Mapping",
                                 cf_transfer_function_type="DIGITAL",
                                 stage_gain=1e6, 
                                 stage_gain_frequency=1.0,
                                 input_units="counts",
                                 output_units="MICROSTRAIN",
                                 input_units_description="Digital counts mapped to int32",
                                 output_units_description="Microstrain units (1e-6 strain)",
                                 numerator=[],
                                 denominator=[],
                                 decimation_factor=1,
                                 decimation_offset=0,
                                 decimation_delay=0,
                                 decimation_correction=0,
                                 decimation_input_sample_rate=ap.get("original_sample_rate")/math.prod(fp["decimation_factor"])
        )
        response_stages.append(microstrain_stage)

    response=Response(
    resource_id="OptoDAS header",
    instrument_sensitivity=InstrumentSensitivity(
        value=ap["sensitivity"],
        frequency=1.0,
        input_units="counts",
        output_units="counts",
        output_units_description="strain per meter",
        input_units_description=ap["sensitivityUnits"]+" phase difference expressed in radians per strain per meter",
        frequency_range_start=0.001,
        frequency_range_end=2000.0,
        frequency_range_db_variation=0,
        ),
        response_stages=response_stages)
    
    sensor=Equipment(
        type="DAS",
        description="OptoDAS",
        manufacturer="ASN",
        serial_number="44",
        installation_date=UTCDateTime(*th.get("net_starttime")),
        resource_id="OptoDAS_44"
        )
    data_logger=Equipment(
        type="DAS",
        description="OptoDAS",
        manufacturer="ASN",
        serial_number="44",
        installation_date=UTCDateTime(*th.get("net_starttime")),
        resource_id="OptoDAS_44"
        )
    pre_amplifier=Equipment(
        type="DAS",
        description="OptoDAS",
        manufacturer="ASN",
        serial_number="44",
        installation_date=UTCDateTime(*th.get("net_starttime")),
        resource_id="OptoDAS_44"
        )
    
    operator=[
        Operator(
            agency="ARDITI",
            contacts=[
                Person(
                    names=["Afonso Loureiro"],
                    agencies=["ARDITI, Funchal, Portugal",
                              "Instituto Dom Luiz, Lisboa, Portugal"],
                    emails=["maloureiro@fc.ul.pt"])],
                    website="www.arditi.pt"
                )
            ]

    inv=Inventory(
        networks=[],
        source=th.get("inventory_source"),
    )

    net=Network(
        code=st[0].stats.network,
        stations=[],
        description=th.get("net_description"),
        start_date=UTCDateTime(ap.get("start_date")),
        end_date=UTCDateTime(ap.get("end_date")),
        total_number_of_stations=th.get("net_no_stations"),
        selected_number_of_stations=len(st),
    )

    for tr in st:
        station_code=tr.stats.station
        if ep.get("encode_station_name"):
            das_channel = tr.stats['decoded_station']
        else:
            das_channel = tr.stats.station
    
        sta=Station(
                code=station_code,
                latitude=tr.stats.coordinates["latitude"],
                longitude=tr.stats.coordinates["longitude"],
                elevation=tr.stats.coordinates["elevation"],
                creation_date=UTCDateTime(*th.get("net_starttime")),
                start_date=UTCDateTime(ap.get("start_date")),
                end_date=UTCDateTime(ap.get("end_date")),
                site=Site(name=f'{th.get("sta_name_prefix")} {tr.stats.channel} channel no.{str(das_channel).zfill(5)}'),
                total_number_of_channels=1,
        )
        sta.operators=operator

        azimuth,dip=get_azimuth_dip(int(tr.stats.distance/ap.get("dx")+1),global_params,roi_channel=False)
        cha=Channel(code="HSF",
            location_code=th.get("location"),
            latitude=tr.stats.coordinates["latitude"],
            longitude=tr.stats.coordinates["longitude"],
            elevation=tr.stats.coordinates["elevation"],
            start_date=UTCDateTime(ap.get("start_date")),
            end_date=UTCDateTime(ap.get("end_date")),
            sample_rate=tr.stats.sampling_rate,
            depth=0,
            azimuth=azimuth,
            dip=dip
        )
        cha.response=response
        cha.equipments=[sensor,data_logger, pre_amplifier]
        cha.sensor=sensor
        cha.data_logger=data_logger
        cha.pre_amplifier=pre_amplifier
        sta.channels.append(cha)
        net.stations.append(sta)

    inv.networks.append(net)

    file=os.path.join(pp.get("MSEED_path"),f'RESP.{th.get("network")}.HSF.xml')
    os.makedirs(pp.get("MSEED_path"), exist_ok=True)
    inv.write(file, format="stationxml", validate=True)

    if ep.get("verbose"):
        print(f'StationXML file written to: {file}')

    if ep.get("verbose", False) == 2:
        print(inv)
        inv[0][-1][0].response.plot(min_freq=0.01,unwrap_phase=True)

    return

def degrade_coordinates(lat, lon, precision=3, error_metres=200):
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

def get_cheby2_coeffs(nyquist, freqmax):
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
        order, wn = signal.cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = signal.cheby2(order, rs, wn, btype="low", analog=0, output="ba")
    return b,a

def zerophase_chebyshev_lowpass_filter(trace, freqmax):
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
        trace.data = signal.filtfilt(b, a, trace.data)


def get_azimuth_dip(channel, global_params,roi_channel=False):
    """
    Get azimuth and dip based on the provided channel

    Parameters:
    - channel (int): The channel number for which coordinates are needed.
    - coords_path (str): The path to the CSV file.
    - coords_DAS (str): The name of the CSV file.

    Returns:
    Tuple (float, float): A tuple containing azimuth and dip.
    """
    ep, pp = global_params['exec_params'], global_params['path_params']

    if ep.get("verbose", False) == 2:
        print(f"Getting azimuth and dip for channel {channel} from {pp['coords_path']}{pp['coords_DAS']}")
        
    df = pd.read_csv(os.path.join(pp['coords_path'],pp['coords_DAS']))

    if roi_channel:
        filtered_df = df.loc[df['ROI_CHANNEL'] == channel]
    else:
        filtered_df = df.loc[df['CHANNEL'] == channel]
    
    if len(filtered_df) == 0:
        azi=0
        dip=0
    else:
        azi=filtered_df['AZIMUTH'].values[0]
        dip=filtered_df['DIP'].values[0]
    
    if azi<0:
        azi+=360

    return azi,dip

def get_data_from_S3(bucket_name, endpoint_url, starttime, endtime, local_path, unsigned=False):
    """
    Download DAS raw files from an S3-compatible endpoint for a single day.

    Parameters:
    - bucket_name: name of the S3 bucket.
    - endpoint_url: S3 endpoint URL.
    - starttime, endtime: datetime.datetime objects defining the interval (must be same day).
    - local_path: local directory where files will be saved (organized by hour).
    - unsigned: if True, use unsigned (anonymous) S3 access.

    Returns:
    - True on success, False on failure or invalid input.
    """
    # Local Imports
    import boto3
    import os
    from botocore.config import Config
    from botocore import UNSIGNED

    # 1. Multi-day Guard: Warn and exit if start and end are on different days
    if starttime.date() != endtime.date():
        print(f"WARNING: starttime ({starttime.date()}) and endtime ({endtime.date()}) "
              "are on different days. This function only processes single-day requests.")
        return False

    # 2. Setup S3 Resource
    if unsigned:
        s3_res = boto3.resource('s3', endpoint_url=endpoint_url, config=Config(signature_version=UNSIGNED))
    else:
        s3_res = boto3.resource('s3', endpoint_url=endpoint_url)

    bucket = s3_res.Bucket(bucket_name)

    try:
        prefix = starttime.strftime('%Y%m%d/dphi/%H')
        dir_structure = starttime.strftime('%Y%m%d/dphi/')
        h5start = starttime.strftime('%Y%m%d/dphi/%H%M%S.h5')
        h5end = endtime.strftime('%Y%m%d/dphi/%H%M%S.h5')
        
        files_to_download = []
        last_seen_key = None
        
        # Check if we should skip the buffer file (Midnight check)
        is_midnight = (starttime.hour == 0 and starttime.minute == 0 and starttime.second == 0)

        for obj in bucket.objects.filter(Prefix=prefix):
            if obj.key > h5end:
                break
            
            if h5start <= obj.key <= h5end:
                # Add buffer file ONLY if not midnight AND it's the first file found
                if not files_to_download and last_seen_key and not is_midnight:
                    files_to_download.append(last_seen_key)
                
                files_to_download.append(obj.key)
            
            last_seen_key = obj.key

        # 3. Execution
        if not files_to_download:
            print("No files found for the specified time range.")
            return True

        target_dir = os.path.join(local_path, dir_structure)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        for key in files_to_download:
            target = os.path.join(target_dir, os.path.basename(key))
            if os.path.exists(target):
                print(f"File {target} already exists in {target_dir}. Skipping download.")
                continue
            print(f"Downloading: {key} to {target}")
            bucket.download_file(key, target)

        print(f"Successfully downloaded {len(files_to_download)} files.")
        return True

    except Exception as e:
        print(f"Error: {e}")
        return False


def baleen_whale_detector(st, global_params, detection_mode='single', channel_idx=None, band_detect=(15,25), env_win=0.2,\
                          thresh=4.5, min_dist=8.0, template_path="./finwhale_100Hz.txt", template_fs=100.0, search_win_sec=1.0, corr_thresh=0.15,\
                          spec_source="raw", spec_fmin=10, spec_fmax=40, spec_nperseg=1042, spec_noverlap=900, spec_nfft=2048, spec_dr=15, cmap="jet"):
    """
    Detector for baleen whale calls in DAS or seismic data.

    This function implements a matched-filter / envelope detector workflow tailored
    for fin whale-like signals. It supports single-channel, auto, or all-channel
    detection modes and returns plotting-ready results and summary statistics.

    Parameters:
    - st: ObsPy Stream object containing traces to analyze.
    - global_params: configuration dictionary (unused by default but kept for API compatibility).
    - detection_mode: 'single', 'auto', or 'all'.
    - other kwargs: detector and spectrogram parameters (see function signature).

    Returns:
    - Varies: prints plots and returns detection summaries depending on mode.
    """
    from scipy.signal import butter, filtfilt, find_peaks, spectrogram, resample_poly
    from obspy.signal.cross_correlation import correlate as obspy_cc
    from datetime import timedelta

    # --- USER SETTINGS
    # ---------------------------------------------------------------------------
    # Detection mode:
    # "single" = run on one channel defined by CHANNEL_IDX
    # "auto"   = choose automatically the strongest channel and run on that one
    # "all"    = run on all channels in the stream
    DETECT_MODE = "single"     # "single", "auto", or "all"
    CHANNEL_IDX = 337          # used only when DETECT_MODE = "single"
    
    # Initial detection parameters
    BAND_DETECT = (15, 25)   # Frequency band used for filtering, amplitude calculation, and peak search
    ENV_WIN = 0.2            # Length of the smoothing window applied to the amplitude
    THRESH = 4.5             # Detection threshold relative to the median filtered amplitude
    MIN_DIST = 8.0           # Minimum spacing between consecutive detections to avoid false positives
    
    # Cross-correlation with a matched filter
    TEMPLATE_PATH = "./finwhale_100Hz.txt"
    TEMPLATE_FS = 100.0     # Original sampling rate of the template. If different from stream, it is resampled
    SEARCH_WIN_SEC = 1.0    # Time window used to search for the best correlation peak
    CORR_THRESH = 0.15      # Minimum correlation value required to keep a detection
    
    # Spectrogram parameters
    SPEC_SOURCE = "raw"     # Trace used for the spectrogram, "raw" or "filtered" data
    SPEC_FMIN = 10          # Frequency interval
    SPEC_FMAX = 40
    SPEC_NPERSEG = 1042     # Number of samples per spectrogram window
    SPEC_NOOVERLAP = 900    # Number of overlapping samples between consecutive windows
    SPEC_NFFT = 2048
    SPEC_DR = 15            # Display dynamic range in dB below the maximum plotted level
    CMAP = "jet"
    
    def bandpass(x, fs, lo, hi, order=4):
        nyq = 0.5 * fs
        lo = max(0.001, float(lo))
        hi = min(float(hi), nyq * 0.99)
        b, a = butter(order, [lo / nyq, hi / nyq], btype="band")
        return filtfilt(b, a, x)
    
    def rms(x):
        x = np.asarray(x, dtype=float)
        return np.sqrt(np.mean(x**2))
    
    def print_stats(name, values, fmt="{:.3f}"):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            print(f"{name} Mean:    n/a")
            print(f"{name} Median:  n/a")
            print(f"{name} Std:     n/a")
            print(f"{name} min:     n/a")
            print(f"{name} max:     n/a")
            return
        print(f"{name} Mean:    " + fmt.format(np.mean(values)))
        print(f"{name} Median:  " + fmt.format(np.median(values)))
        print(f"{name} Std:     " + fmt.format(np.std(values)))
        print(f"{name} min:     " + fmt.format(np.min(values)))
        print(f"{name} max:     " + fmt.format(np.max(values)))

    def centers_to_edges(c, start=None, end=None):
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
    
    def trace_label(tr, idx):
        station = getattr(tr.stats, "station", "")
        if station:
            return station
        return f"channel_{idx}"
    
    def prepare_template(fs):
        tpl = np.loadtxt(TEMPLATE_PATH).astype(np.float64)
        tpl -= tpl.mean()
        if TEMPLATE_FS != fs:
            tpl = resample_poly(
                tpl,
                up=int(round(fs)),
                down=int(round(TEMPLATE_FS))
            ).astype(np.float64)
            tpl -= tpl.mean()
        n_tpl = len(tpl)
        pulse_width = (n_tpl - 1) / fs
        half_pulse_samp = int(round(0.5 * pulse_width * fs))
        n_search_window = int(round(SEARCH_WIN_SEC * fs))
        # zero + master + zero
        tpl_plus = np.concatenate([np.zeros(n_tpl), tpl, np.zeros(n_tpl)])
        return tpl, tpl_plus, n_tpl, half_pulse_samp, n_search_window
    
    def auto_choose_channel(stream):
        peaks = []
        for i, tr in enumerate(stream):
            x_raw = tr.data.astype(np.float64)
            fs = float(tr.stats.sampling_rate)
            x_det = bandpass(x_raw, fs, *BAND_DETECT)
            peaks.append(np.nanmax(np.abs(x_det)))
        return int(np.nanargmax(peaks))

    def process_trace(tr, idx):
        x_raw = tr.data.astype(np.float64)
        fs = float(tr.stats.sampling_rate)
        start_time = tr.stats.starttime.datetime
        channel_label = trace_label(tr, idx)
        tpl, tpl_plus, n_tpl, half_pulse_samp, n_search_window = prepare_template(fs)
        # Detection trace
        x_det = bandpass(x_raw, fs, *BAND_DETECT)
        env = np.abs(x_det)
        win_n = max(1, int(ENV_WIN * fs))
        amp_smooth = np.convolve(env, np.ones(win_n) / win_n, mode="same")
        thr = THRESH * np.median(env)
        dist = max(1, int(MIN_DIST * fs))
        raw_picks, _ = find_peaks(
            amp_smooth,
            height=thr,
            distance=dist
        )
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
            
            cc = obspy_cc(
                st2,
                tpl_plus,
                n_search_window,
                normalize="naive",
                method="direct",
                demean=True
            )
            i_max = int(np.argmax(cc))
            master_corr = float(np.max(cc))
            master_shift = (i_max - n_search_window) / fs
            if master_corr < CORR_THRESH:
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
            signal_rms = float(rms(signal_win))
            noise_rms = float(rms(noise_win))
            snr = np.nan if noise_rms == 0 else signal_rms / noise_rms
            tsec = refined_start / fs
            abs_time = start_time + timedelta(seconds=float(tsec))
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
            "idx": idx,
            "label": channel_label,
            "fs": fs,
            "start_time": start_time,
            "x_raw": x_raw,
            "x_det": x_det,
            "amp_smooth": amp_smooth,
            "thr": thr,
            "raw_picks": raw_picks,
            "results": results,
        }
    def choose_plot_result(processed):
        if len(processed) == 1:
            return processed[0]
        # Prefer the channel with the most kept detections
        kept_counts = [len(p["results"]) for p in processed]
        if np.max(kept_counts) > 0:
            return processed[int(np.argmax(kept_counts))]
        # If none kept, use the one with most raw picks
        raw_counts = [len(p["raw_picks"]) for p in processed]
        return processed[int(np.argmax(raw_counts))]

    # --- BASIC CHECKS
    # ---------------------------------------------------------------------------
    if st==None or len(st) == 0:
        raise RuntimeError("Stream 'st' is empty.")
    
    if not os.path.exists(TEMPLATE_PATH):
        raise FileNotFoundError(f"Template not found: {TEMPLATE_PATH}")
    
    # --- CHANNEL SELECTION
    # ---------------------------------------------------------------------------
    if DETECT_MODE not in ("single", "auto", "all"):
        raise ValueError("DETECT_MODE must be 'single', 'auto', or 'all'")
    
    if DETECT_MODE == "single":
        if CHANNEL_IDX is None:
            raise ValueError("Set CHANNEL_IDX when DETECT_MODE = 'single'")
        if CHANNEL_IDX < 0 or CHANNEL_IDX >= len(st):
            raise IndexError(f"CHANNEL_IDX out of range: {CHANNEL_IDX}")
        process_indices = [CHANNEL_IDX]
    
    elif DETECT_MODE == "auto":
        auto_idx = auto_choose_channel(st)
        process_indices = [auto_idx]
    
    else:  # DETECT_MODE == "all"
        process_indices = list(range(len(st)))
    
    
    # --- RUN DETECTOR
    # ---------------------------------------------------------------------------
    processed = []
    for idx in process_indices:
        processed.append(process_trace(st[idx], idx))
    
    # Channel used for plotting
    plot_result = choose_plot_result(processed)
    
    # Flatten all kept detections for global summary
    all_results = []
    for p in processed:
        for rec in p["results"]:
            rec2 = rec.copy()
            rec2["channel_idx"] = p["idx"]
            rec2["channel_label"] = p["label"]
            all_results.append(rec2)
    
    
    # --- OUTPUT
    # ---------------------------------------------------------------------------
    if DETECT_MODE in ("single", "auto"):
        p = processed[0]
        print(f"Detection mode       : {DETECT_MODE}")
        print(f"Chosen channel index : {p['idx']}")
        print(f"Chosen channel label : {p['label']}")
        print(f"Sampling rate        : {p['fs']:.1f} Hz")
        print(f"Frequency BW         : {BAND_DETECT} Hz")
        print(f"ENV_WIN              : {ENV_WIN} s")
        print(f"Raw picks            : {len(p['raw_picks'])}")
        print(f"Kept detections      : {len(p['results'])}")
    
        for i, rec in enumerate(p["results"], start=1):
            print(
                f"{i:02d} | t = {rec['time_s']:7.2f} s | {rec['datetime']} | "
                f"Cmax = {rec['corr']:.3f} | Shift = {rec['shift_s']:.3f} s | "
                f"Amp = {rec['amp_max']:.3f} | Signal RMS = {rec['signal_rms']:.3f} | "
                f"Noise RMS = {rec['noise_rms']:.3f} | SNR = {rec['snr']:.2f}"
            )
    
    else:
        print(f"Detection mode       : {DETECT_MODE}")
        print(f"Processed channels   : {len(processed)}")
        print(f"Frequency BW         : {BAND_DETECT} Hz")
        print(f"ENV_WIN              : {ENV_WIN} s")
        print()
    
        nonzero = [p for p in processed if (len(p["raw_picks"]) > 0 or len(p["results"]) > 0)]
        print(f"Channels with picks  : {len(nonzero)}")
        print(f"Total kept detections: {len(all_results)}")
        print()
    
        # Print concise summary only for channels with at least one raw or kept detection
        for p in nonzero:
            print(
                f"idx={p['idx']:5d} | label={p['label']} | "
                f"raw={len(p['raw_picks']):3d} | kept={len(p['results']):3d}"
            )
    
    # --- SUMMARY STATS
    print("*******************************")
    print(f"Number of lines found {len(all_results)}")
    
    if len(all_results) > 0:
        cmax_arr = [r["corr"] for r in all_results]
        shift_arr = [r["shift_s"] for r in all_results]
        snr_arr = [r["snr"] for r in all_results]
    
        print_stats("Cmax", cmax_arr, "{:.3f}")
        print()
        print_stats("Shift", shift_arr, "{:.3f}")
        print()
        print_stats("SNR", snr_arr, "{:.3f}")
    
    
    # --- PLOT: DETECTION RESULTS + SPECTROGRAM
    # ---------------------------------------------------------------------------
    picks = np.array([rec["sample"] for rec in plot_result["results"]], dtype=int)
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
    x_spec = x_raw if SPEC_SOURCE == "raw" else x_det
    
    f, tt, Sxx = spectrogram(
        x_spec,
        fs=fs,
        window="hann",
        nperseg=SPEC_NPERSEG,
        noverlap=SPEC_NOOVERLAP,
        nfft=SPEC_NFFT,
        detrend=False,
        scaling="spectrum",
        mode="magnitude"
    )
    
    Sxx_db = 20 * np.log10(Sxx + 1e-12)
    sel = (f >= SPEC_FMIN) & (f <= SPEC_FMAX)
    
    duration = len(x_spec) / fs
    t_edges = centers_to_edges(tt, start=0.0, end=duration)
    f_sel = f[sel]
    f_edges = centers_to_edges(f_sel)
    
    vmax = np.percentile(Sxx_db[sel, :], 99.9)
    vmin = vmax - SPEC_DR
    
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(12, 6), sharex=True,
        gridspec_kw={"height_ratios": [1, 2]}
    )
    fig.subplots_adjust(right=0.80, hspace=0.10)
    
    # Top panel: smoothed amplitude
    ax0.plot(t, amp_smooth, lw=1.0, color="black", label="Smoothed amplitude")
    ax0.axhline(thr, color="red", ls="--", lw=1.0, label="Amplitude threshold")
    
    if len(raw_picks):
        ax0.plot(
            raw_picks / fs,
            amp_smooth[raw_picks],
            "o",
            ms=4,
            mfc="none",
            mec="orange",
            label="Raw picks"
        )
    
    ax0.set_ylabel("Amplitude")
    ax0.set_title(f"Channel {channel_idx} ({channel_label}) | start = {start_time}")
    ax0.set_xlim(0, duration)
    ax0.legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        frameon=False
    )
    
    # Bottom panel: spectrogram
    ax1.pcolormesh(
        t_edges,
        f_edges,
        Sxx_db[sel, :],
        shading="flat",
        cmap=CMAP,
        vmin=vmin,
        vmax=vmax
    )
    
    ax1.set_ylim(SPEC_FMIN, SPEC_FMAX)
    ax1.set_xlim(0, duration)
    ax1.set_ylabel("Frequency (Hz)")
    ax1.set_xlabel("Time (s)")
    
    for p in picks:
        ax1.axvline(p / fs, color="white", ls="--", lw=1.0, alpha=0.9)
    
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.savefig('Whale_detector.png', bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close()
