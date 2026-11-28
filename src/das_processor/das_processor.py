"""DAS processing pipeline: the ``DASProcessor`` class and ``process_events``.

Supporting code lives in sibling modules (``das_io``, ``coords``, ``stationxml``, ``filters``,
``scaling``, ``s3``, ``whale``, ``plotting``, ``_log``). Their public names are re-exported here, so
``from .das_processor import get_coordinates`` and similar imports keep working.
"""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import datetime
import inspect
import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

# --- Third party -------------------------------------------------------------
import h5py
import numpy as np
import pandas as pd
import scipy.fft as sfft
import scipy.signal as sps
from obspy import Stream, Trace, UTCDateTime
from obspy.core.util import AttribDict
from scipy import ndimage
from tqdm.auto import tqdm

# --- Local -------------------------------------------------------------------
# Names moved to sibling modules are re-exported here for backwards compatibility.
from ._log import _setup_logging, TRACE
from .coords import (  # noqa: F401
    clear_coords_cache, degrade_coordinates, encode_station_name, get_azimuth_dip,
    get_coordinates, get_travel_times,
)
from .das_io import (  # noqa: F401
    convert_index_to_rel_time, create_time_axis, DASDataFrame, DictFile, find_DAS_files,
    get_data_indexes, get_filemeta, load, load_DAS_files, load_keys, save, save_to_DAS_file,
    str2datetime, unwrap, wrap,
)
from .filters import (  # noqa: F401
    design_decimation_fir, get_cheby2_coeffs, zerophase_chebyshev_lowpass_filter,
    zerophase_effective_fir,
)
from .plotting import PlottingMixin
from .s3 import get_data_from_S3  # noqa: F401
from .scaling import check_int32_scaling, INT32_MAX, INT32_SCALES  # noqa: F401
from .stationxml import station_xml_creation  # noqa: F401
from .whale import detect_baleen_whales, DetectionParams, SpectrogramParams

logger = logging.getLogger(__name__)


# --- Standard library --------------------------------------------------------

# --- Third party -------------------------------------------------------------

# --- Local -------------------------------------------------------------------
# Names moved to sibling modules are re-exported here for backwards compatibility.



###########################################
############# Helper functions ############
###########################################


def create_global_params(defaults: Any, sources: Any, exec_params: dict[str, Any]) -> dict[str, Any]:
    """Constructs and returns the global configuration dictionary."""
    return {
        "freq_params": defaults.freq_params,
        "exec_params": exec_params,
        "path_params": defaults.path_params,
        "trace_header": defaults.trace_header,
        "channel_ranges": defaults.channel_ranges,
        "offset_ranges": defaults.offset_ranges,
        "sources": sources.sources,
        "station_coordinates": defaults.station_coordinates,
        "acquisition_params": defaults.acquisition_params,
        "clock_offset": defaults.clock_offset,
    }

def validate_sources(sources) -> dict:
    """
    Validates a dictionary of source configurations against defined rules:
    1. Must have 'start_time' and either 'window_duration' or 'end_time'.
    2. Entries with type 'seismic' must have 'magnitude', 'latitude', 'longitude', and 'depth'.
    3. All entries except 's3_download' and 'delete_data' mode must have 'channel_ranges'.
    4. All entries must have 'type' and 'mode'.
    
    Returns a dictionary of error messages for each entry.
    """
    validation_errors = {}

    for name, config in sources.sources.items():
        errors = []
        
        # Requirement 1: Each entry must have 'start_time' and either 'window_duration' or 'end_time'
        if "start_time" not in config:
            errors.append("Missing required field: 'start_time'")
        
        has_window = "window_duration" in config
        has_end = "end_time" in config
        if not (has_window or has_end):
            errors.append("Must have either 'window_duration' or 'end_time'")
            
        # Requirement 2: Entries of type 'seismic' must have magnitude, latitude, longitude, and depth
        if config.get("type") == "seismic":
            for field in ["magnitude", "latitude", "longitude", "depth"]:
                if field not in config:
                    errors.append(f"Seismic entry missing required field: '{field}'")
                    
        # Requirement 3: All entries except 's3_download' and 'delete_data' must have 'channel_ranges' and 'type'
        if config.get("mode") not in ("s3_download", "delete_data"):
            if "type" not in config:
                errors.append("Missing required field: 'type'")
            if "channel_ranges" not in config:
                errors.append("Missing required field: 'channel_ranges'")
        
        # Requirement 4: All entries must have 'mode'
        if "mode" not in config:
            errors.append("Missing required field: 'mode'")
                
        if len(errors) != 0:
            validation_errors[name] = errors

    return validation_errors

def validate_source_channel_ranges(
    events_to_process: list[str], global_params: dict
) -> None:
    excluded_modes = {"time_lapse", "s3_download", "delete_data"}
    sources = global_params.get("sources", {})
    channel_ranges = global_params.get("channel_ranges", {})

    for event in events_to_process:
        if event not in sources:
            continue
        source = sources[event]
        event_mode = source.get("mode")
        if event_mode in excluded_modes:
            continue

        source_type = source.get("channel_ranges")
        if source_type not in channel_ranges:
            raise ValueError(
                f"{event}: '{source_type}' does not exist in the defined channel_ranges."
            )

def _validate_seismic_event_fields(
    events_to_process: list[str], global_params: dict
) -> None:
    """Require seismic source metadata only for selected seismic events."""
    required_fields = ("magnitude", "latitude", "longitude", "depth")
    sources = global_params.get("sources", {})

    for event in events_to_process:
        source = sources.get(event)
        if source is None or source.get("type") != "seismic":
            continue

        missing_fields = [
            field
            for field in required_fields
            if source.get(field) is None
            or (isinstance(source.get(field), str) and not source[field].strip())
        ]
        if missing_fields:
            raise ValueError(
                f"{event}: seismic events require "
                f"{', '.join(missing_fields)}."
            )


        


########################################
###### Event processing functions ######
########################################

def process_events(events_to_process: list[str], global_params: dict[str, Any]):
    _setup_logging(verbose=global_params.get("exec_params", {}).get("verbose"))
    
    sources = global_params.get("sources", {})
    station = global_params.get("exec_params", {}).get("station", "")   
    station_coordinates = global_params.get("station_coordinates", {})
    exec_params = global_params.get("exec_params", {})

    # Execution parameters go to the log file, and to the console when verbose is enabled
    logger.debug("EXEC Parameters are:")
    for key, value in exec_params.items():
        if value:
            logger.debug(f"{key}: {value}")
    logger.debug('#############################')

    validate_source_channel_ranges(events_to_process, global_params)
    _validate_seismic_event_fields(events_to_process, global_params)

    if not any(event in sources for event in events_to_process):
        logger.error("List of events to process did not match any event in sources.py.")
        logger.info("Nothing to do in `process_events`")
        logger.warning(f"List of events to process ({events_to_process}) did not match any event in sources.py.")
        logger.warning("Nothing to process.")
        return None
    
    logger.info("######################")
    logger.info('»»» Processing event windows.')
    for event in tqdm(events_to_process, desc="Processing events"):
        if not event in sources:
            valid_events = "\n".join(key for key in sources)
            logger.warning(f'{event} not in list of defined events.\nValid events are:\n{valid_events}')
            continue
        else:
            logger.info(f'Working on event: {event}.')

        logger.info(f'#########################\n#########################\n\n')
        logger.info(f'Processing event {event}. station: {station},station_coordinates: {station_coordinates[station]}')

        start_time_event = time.time()

        # 1. Initialize
        DAS = DASProcessor(global_params)
        DAS.event_info = sources.get(event, {})

        # Get the operation mode, defaulting to 'seismic' if not specified
        event_mode = DAS.event_info.get("mode", "")
        logger.info(f"Event mode: {event_mode}")
        success = False

        with DAS.temporary_event_paths(event):
            if event_mode == "time_lapse":
                logger.info(f"Processing time lapse for event: {event}")
                DAS.process_time_lapse(event)

            elif event_mode == "s3_download":
                logger.info(f"Downloading data from S3 for event: {event}")
                success=get_data_from_S3(
                bucket_name=DAS.event_info["s3_bucket"],
                endpoint_url=DAS.event_info["endpoint_url"],
                starttime=UTCDateTime(DAS.event_info.get("start_time")),
                endtime=UTCDateTime(DAS.event_info.get("end_time")),
                local_path=DAS.event_info["DAS_data_path"],
                unsigned=DAS.event_info["unsigned", False]
            )
                if success:
                    logger.info(f"S3 download for {event} completed successfully.")
                    success = False
            elif event_mode == "delete_data":
                logger.info(f"Deleting data for event: {event}")
                DAS.delete_event_data(event)
            else:
                # 2. Load event (check if files exist, load data, populate metadata)
                logger.info(f"Loading event data for event: {event} in mode: {event_mode}")
                success = DAS.load_event_data(event, mode=event_mode)

            # 3. Check result
            if success:
                logger.info("Event data loaded successfully!"
                      f"Shape of DAS data (samples x channels): {DAS.dfdas.shape}"
                      f"DFDAS start time: {DAS.datetime_start}, duration: {DAS.duration} seconds.")
                logger.info(f"Event data loaded successfully for event: {event}\n"
                             f"Shape of DAS data (samples x channels): {DAS.dfdas.shape}\n"
                             f"DFDAS start time: {DAS.datetime_start}, duration: {DAS.duration} seconds.")

                event_info = DAS.sources.get(event, {})
                process_chain = event_info.get('process_chain')
                if process_chain: # This evaluates to True if it's a non-empty list
                    logger.info(f"Executing process chain for {event}...")
                    logger.info(f"Executing process chain for {event}...")
                    DAS.run_process_chain() 
                    logger.info("Process chain completed.")
                else:
                    logger.info(f"No process chain defined for {event}. Skipping automated steps.")
                    logger.warning(f"No process chain defined for {event}. Skipping automated steps.")
            elif event_mode != "time_lapse":
                logger.error(f"Event processing failed or no data was found for event: {event}")
                logger.error("Event processing failed or no data was found.")

            end_time_event = time.time()
            elapsed_time = end_time_event - start_time_event

        logger.info(f"\n>>> FINISHED: {event}")
        logger.info(f">>> Execution time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        logger.info("#########################\n")

        if event_mode == "time_lapse" or event_mode == "s3_download" or event_mode == "delete_data":
            logger.info(f'Processing mode "{event_mode}" does not require data loading.\nNo data is available for further processing.')
        else:
            logger.info(f"Loaded data for {event} was returned.")

    if "DAS" not in locals() or DAS is None:
        logger.warning("No data loaded.")
        logger.warning("No data loaded.")
    return DAS

# define class for storing results

class DASProcessor(PlottingMixin):
    def __init__(self, global_params: dict[str, Any]) -> None:
        self.exec_params = global_params.get('exec_params', {})
        self.path_params = global_params.get('path_params', {})
        self.sources = global_params.get('sources', {})
        self.channel_ranges = global_params.get('channel_ranges', {})
        self.offset_ranges = global_params.get('offset_ranges', {})
        self.station_coordinates = global_params.get('station_coordinates', {})
        self.acquisition_params = global_params.get('acquisition_params', {})
        self.trace_header = global_params.get('trace_header', {})
        self.freq_params = global_params.get('freq_params', {})
        self.clock_offset = global_params.get('clock_offset',0)

        self.stationxml_params = {}
        self.dfdas_decimation_steps = []    
        self.network = self.exec_params.get("network", "XX")  # Default to '5V' if not specified

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

        self.whale_figure = None   # last figure drawn by baleen_whale_detector

        # Set up basic logging configuration with execution timestamp
        _setup_logging(verbose=self.exec_params.get("verbose"))


    @staticmethod
    def _make_step(
        factor: int, input_rate: float, *,
        b: np.ndarray | None = None, a: np.ndarray | None = None,
        fir: np.ndarray | None = None, label: str | None = None,
    ) -> dict[str, Any]:
        if (fir is None) == (b is None or a is None):
            raise ValueError("Provide either fir, or both b and a.")
        return {
            "factor": int(factor),
            "input_sampling_rate": float(input_rate),
            "fir": None if fir is None else np.asarray(fir),
            "b": None if b is None else np.asarray(b),
            "a": None if a is None else np.asarray(a),
            "label": label,
        }

    def get_decimation_steps(self) -> list:
        """Steps applied to self.dfdas, then those applied to the current stream."""
        st_steps = (getattr(self.mseed_stream, "decimation_steps", [])
        if self.mseed_stream is not None else [])
        return list(self.dfdas_decimation_steps) + list(st_steps)
    
    def get_current_params(self) -> dict[str, Any]:
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
            'stationxml_params': {**self.stationxml_params,
                                  "decimation_steps": self.get_decimation_steps()},
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
                    logger.info(f"Override: Setting path_params[{key}] to {new_val}")
                    self.path_params[key] = new_val

            # Apply execution overrides (using 'is not None' so False is successfully set!)
            for key, new_val in exec_overrides.items():
                if new_val is not None:
                    logger.info(f"Override: Setting exec_params[{key}] to {new_val}")
                    self.exec_params[key] = new_val

            yield  # <--- Processing runs here

        finally:
            # 4. RESET to original defaults securely
            for key, old_val in original_paths.items():
                self.path_params[key] = old_val

            for key, old_val in original_exec.items():
                self.exec_params[key] = old_val

            logger.info("Paths and execution parameters restored to defaults.")

    def run_process_chain(self, suffix=""):
        """
            Iterates through the 'process_chain' defined in the source
            and executes methods dynamically.
            Passes a suffix (e.g., '_timelapse_1') to methods for unique filenames.
            The stored chain definition is never modified, so this method can be
            called repeatedly with different suffixes.
        """

        chain = self.event_info.get('process_chain', [])

        for step in chain:
            operation = step.get('operation')

            # Work on a copy: the original step['params'] stays untouched
            params = dict(step.get('params') or {})

            logger.info(f"Processing step: {operation} with params: {params}")
            
            if not hasattr(self, operation):
                logger.warning(f"Warning: Method {operation} not found in DASProcessor.")
                logger.warning(f"Method {operation} not found in DASProcessor.")
                continue

            method = getattr(self, operation)
            logger.info(f"Executing step: {operation} with params: {params}")
            logger.debug(f"--- Executing: {operation} with params {params} ---")
            
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
            logger.info(f"Completed step: {operation}")

    def get_time_lapse_windows(self):
        """Generates start times for continuous windows."""
        start_time = UTCDateTime(self.event_info.get("start_time"))
        end_time = UTCDateTime(self.event_info.get("end_time"))
        duration = self.event_info.get("window_duration")

        start_times = []
        current = start_time
        while current < end_time:
            start_times.append(current)
            current += datetime.timedelta(seconds=duration)

        logger.info(f'{len(start_times)} windows of {duration}s to process')
        return start_times, duration
    
    def process_time_lapse(self, event):
        """Main loop for time-lapse processing."""
        logger.info(f'»»» Working on time lapse: {event}')
        
        # 1. Generate the windows
        start_times, window_duration = self.get_time_lapse_windows()
        self.duration = window_duration # Set window length for the processor

        events_path = self.path_params.get("events_path")
        event_safe = event.replace(" ", "_")

        # 2. Iterate through windows
        for i, start_dt in enumerate(tqdm(start_times, desc="Time lapse windows")):
            suffix = f"_{start_dt.strftime('%Y%m%d_%H%M')}"
            target_hdf5 = os.path.join(events_path, f"{event_safe}{suffix}.hdf5")

            # Check if we should skip before loading raw data
            if os.path.exists(target_hdf5):
                logger.info(f"Window {i+1}/{len(start_times)}: Skipping (Already exists: {os.path.basename(target_hdf5)})")
                continue

            logger.info(f"\n--- Processing Window {i+1}/{len(start_times)}: {start_dt} ---")
            
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
                logger.warning(f"Skipping window {i+1} due to loading error.")
    
        return False

    def delete_event_data(self, event, dry_run=True):
        """Delete raw DAS files overlapping the event's configured time range."""

        start_time=UTCDateTime(self.event_info.get("start_time")),

        if self.event_info.get("end_time") is None and self.event_info.get("duration") is not None:
            duration_event = self.event_info.get("duration")
            if not isinstance(duration_event, (int, float)):
                raise ValueError(f"Event {event!r} has an invalid 'duration' value.")
            end_time = start_time + datetime.timedelta(seconds=duration_event)    
        else:
            end_time = UTCDateTime(self.event_info.get("end_time"))
        
        if start_time is None or end_time is None:
            raise ValueError(
                f"Event {event!r} must define both 'start_time' and 'end_time' or 'start_time' and 'duration' "
                "to delete its DAS data."
            )

        duration = (end_time - start_time).total_seconds()
        if duration <= 0:
            raise ValueError(f"Event {event!r} has an invalid time range.")

        experiment_path = self.path_params.get("experiment_path")
        if not experiment_path:
            raise ValueError("The DAS experiment path is not configured.")

        files, _, _ = find_DAS_files(
            experiment_path,
            start=start_time,
            duration=duration,
            show_header_info=self.exec_params.get("verbose", False),
            load_file_from_start=False,
        )

        warning = f"WARNING: deleting raw DAS data for event {event}."
        logger.warning("%s", warning)
        logger.warning(warning)
        logger.info(f"Deletion summary:\n"
            f"  Time range: {start_time} to {end_time}\n"
            f"  Experiment path: {experiment_path}\n"
            f"  {len(files)} files selected :")
        for filepath in files:
            logger.info(f"  - {filepath}")

        if not files:
            logger.warning("No matching DAS files found; nothing was deleted.")
            logger.info(f"No DAS files found for event {event}; nothing was deleted.")
            return []

        if dry_run:
            logger.info("Dry run enabled; no files were deleted.")
            logger.info(f"Dry run: no DAS files deleted for event {event}.")
        else:
            for filepath in files:
                os.remove(filepath)
            message=f"Deleted {len(files)} DAS files for event {event}."
            logger.info("%s", message)
            logger.info(message)

        return

    def load_event_data(self, event, mode='seismic', replace_DAS=False):
        """
        The main pipeline execution method. 
        """
        logger.info(f'Working on event: {event}.')

        # 1. Initial setup
        self.dfdas = None
        self.dfdas_offsets = None
        self.mseed_stream = None
        self.event = event
        station = self.exec_params.get("station")
        
        # reset decimation history
        self._reset_decimation_history(domain='dfdas')

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
            logger.error(f"No files found for {event}. Execution aborted for this event.")
            return False 
        
        # 4. Handle Stale GPS Logic
        # We check if the GPS status changed after looking at the file headers
        if self.acquisition_params.get("gps_status") == True:
            if self.sources.get(event, {}).get("GPS") == False:
                logger.warning("GPS information was stale. Recalculating start times.")
                figsize, self.datetime_start, self.duration = self.calculate_trim_times(event, station)
                # Re-locate files with updated time
                self.locate_event_data(event)

        # 5. Check if we found anything
        if not self.files:
            logger.warning(f'No data available for event {event}. Skipping!')
            return False # Signal to the loop to continue

        # 6. Calculate Channels
        self.channels, min_ch, max_ch, ch_step = self.calculate_channels(event)

        # Exit early if channels array is empty
        if self.channels.size == 0:  # Use .size for NumPy arrays
            logger.error(f"No valid channels found for event {event}. Exiting.")
            return False

        # 7. Final Load
        self.load_dfdas_file(event, replace_DAS=replace_DAS)

        logger.debug(f"LOAD_EVENT_DATA: Event {event} loaded successfully with {len(self.channels)} channels.")
        
        return True # Success!

    def locate_event_data(self,event):

        events_path = self.path_params.get("events_path")
        experiment_path = self.path_params.get("experiment_path")
        # replace spaces with underscores in event name for filename
        hdf5_file = f'{events_path}/{event.replace(" ", "_")}.hdf5'

        self.event=event

        logger.info(f'Finding files for event: {event}')
        logger.info(f'Trying file list from {events_path}')

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
            logger.info(f'Local HDF5 not found. Searching raw data in {experiment_path}')
            date_start= self.datetime_start.strftime("%Y%m%d")
            target_path=Path(experiment_path) / date_start
            if not target_path.exists():
                logger.warning(f"Skipping event {event} due to missing directories.")
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
        files, chIndex, samples = find_DAS_files(
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
                ap["gauge_length"] = hf_in['/header/gaugeLength'][()]
                ap["dx"] = hf_in['/header/dx'][()]
                ap["original_sample_rate"] = 1.0 / hf_in['/header/dt'][()]
                ap["sensitivity"] = hf_in['/header/sensitivities'][0][0]
                ap["sensitivityUnits"] = hf_in['/header/sensitivityUnits'][0].decode('utf-8')
                
                # Calculate overlap
                ap["overlap"] = ap["roi_step"] * ap["dx"] / ap["gauge_length"]

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

        if einfo.get("window_duration") is None:
            end_time = einfo.get("end_time")
            start_time = einfo.get("start_time")
            einfo["window_duration"] = (end_time - start_time).total_seconds()
            logger.warning(f"`window_duration` was not set for {event}. Calculated from `start_time` and `end_time`.")
            logger.debug(f"`window_duration` was not set for {event}. Calculated from `start_time` and `end_time`.")

        if self.roi_step is None:
            logger.warning("ROI step is not set. Attempting to retrieve from loaded data.")
            self.roi_step = self.filemeta['demodSpec']['roiDec'][0]
            if self.roi_step is None:
                logger.warning("ROI step could not be retrieved from current dfdas data."
                                 " Attempting to retrieve from files.")
                experiment_path = self.path_params.get("experiment_path")
                start_time = UTCDateTime(einfo.get("start_time")).datetime
                self.get_roi_step(experiment_path, start_time)
            if self.roi_step is None:
                logger.error("ROI step could not be determined. Please check the data files.")
                raise ValueError("ROI step is not set and could not be retrieved from files.")

        logger.info(f'Processing event {event}, of the type {einfo.get("type")}')
        logger.info(f'Processing event {event}, of the type {einfo.get("type")}')
        logger.info('Calculating trim times and figure sizes...')
        logger.info('Calculating trim times and figure sizes.')

        # Determine clock offset based on GPS status stored in the class
        if ap.get("gps_status"):
            self.clock_offset = 0
        else:
            if self.clock_offset == 0:
                logger.warning("No GPS timing reported, but clock was not corrected.")
                logger.warning(f"No GPS timing reported, but clock was not corrected.\n"
                      f"Revise 'clock_offset' parameter in defaults.py")

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
            trim_time = base_time + arrivals_p[0].time - einfo.get("window_in") + datetime.timedelta(seconds=self.clock_offset)
            logger.info(f'Time offset correction of {self.clock_offset}s applied.')
        else:
            self.figsize = (10, 8)
            trim_time = UTCDateTime(einfo.get("start_time")) - einfo.get("window_in")
            logger.info("Time offset not corrected.")

        # Store results in the class instance
        self.datetime_start = trim_time.datetime
        self.duration = datetime.timedelta(seconds=einfo.get("window_duration") + einfo.get("window_in"))

        logger.info(f'Datetime start: {self.datetime_start}')
        return self.figsize, self.datetime_start, self.duration

    def trim_DAS_data(self, event):
        """
        Trims the DAS data to the specified time window.
        """
        window_in = pd.to_timedelta(self.event_info.get("window_in"), unit='s')
        window_duration = pd.to_timedelta(self.event_info.get("window_duration"), unit='s')
        start_time = pd.to_datetime(self.event_info.get("start_time"))

        self.dfdas = self.dfdas.truncate(before=start_time - window_in, after=start_time + window_duration)
        logger.info("%s %s %s %s", 'DAS data trimmed from ', start_time - window_in, ' to ', start_time + window_duration)
        logger.info("%s %s", 'with a start time of ', start_time)

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
            self.files, ch_idx, samples = find_DAS_files(
                pp.get("experiment_path"), channels=channels, 
                start=datetime_start, duration=duration, load_file_from_start=False
            )
            self.dfdas = load_DAS_files(self.files, ch_idx, samples, integrate=False)
        else:
            if os.path.exists(local_hdf5):
                try:
                    self.files = [local_hdf5]
                    # This avoids channel re-numbering
                    ch_idx, samples = get_data_indexes(self.files[0])
                    logger.debug("%s", self.files[:])
                    self.dfdas = load_DAS_files(self.files, ch_idx, samples, integrate=False)
                except ValueError:
                    # Fallback to experiment path if local HDF5 fails
                    self.files, ch_idx, samples = find_DAS_files(
                        pp.get("experiment_path"), channels=channels, 
                        start=datetime_start, duration=duration, load_file_from_start=False
                    )
                    self.dfdas = load_DAS_files(self.files, ch_idx, samples, integrate=False)
            else:
                self.files, ch_idx, samples = find_DAS_files(
                    pp.get("experiment_path"), channels=channels, 
                    start=datetime_start, duration=duration, load_file_from_start=False
                )
                self.dfdas = load_DAS_files(self.files, ch_idx, samples, integrate=False)

        self.dfdas = DASDataFrame(self.dfdas, meta=self.dfdas.meta)

        # 2. Sanity check using acquisition_params (updated by previous methods)
        if einfo.get("GPS") and ap.get('gps_status') != 0:
            logger.warning("Event defined as having locked GPS, but file header indicates otherwise.\nClock offset should be calculated.")

        logger.debug(f'chIndex: {ch_idx}')
        logger.debug(f'columns: {self.dfdas.columns}')

        # 3. Store Metadata
        self.filemeta = get_filemeta(self.files[0], 2)

        # 3.1 update self.roi_step in case it was not set before
        if self.roi_step is None:
            logger.warning("ROI step was not set upon loading. Retrieving from file metadata.")
            self.roi_step = self.filemeta['demodSpec']['roiDec'][0]
            ap["roi_step"] = self.roi_step

        # 3.2 Update self.first_chan and self.last_chan if they were not set before
        if self.first_chan is None:
            logger.warning("First channel was not set upon loading. Retrieving from file metadata.")
            self.first_chan = self.filemeta['header']['channels'][0]
            ap["first_chan"] = self.first_chan
        if self.last_chan is None:
            logger.warning("Last channel was not set upon loading. Retrieving from file metadata.")
            self.last_chan = self.filemeta['header']['channels'][-1]
            ap["last_chan"] = self.last_chan

        # 4. Save if required
        if ep.get("save_DFDAS") and event_mode != "time_lapse":
            os.makedirs(pp.get("events_path"), exist_ok=True)
            file_exists = os.path.exists(local_hdf5)
            
            if file_exists and not replace_DAS:
                logger.info(f'Local HDF5 already exists for event {event}. Skipping save.')
            else:
                if file_exists:
                    logger.info(f'Local HDF5 already exists for event {event}. Overwriting.')

                filename_out = save_to_DAS_file(
                    self.dfdas, local_hdf5, "processed", **self.filemeta
                    )
                logger.info(f'Window for event {event} saved to {filename_out[0]}')

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
            logger.warning("No 't' column found in DataFrame for gap checking.")
            return []

        time_diffs = self.dfdas['t'].diff().dropna()
        expected_dt = self.acquisition_params.get("roi_step", 1) * self.acquisition_params.get("dx", 1) / self.acquisition_params.get("original_sample_rate", 1)
        
        # Define a threshold for identifying gaps (e.g., 1.5 times the expected dt)
        gap_threshold = expected_dt * 1.5
        gaps = time_diffs[time_diffs > gap_threshold]

        if not gaps.empty:
            logger.warning(f"Warning: Detected {len(gaps)} gaps in the data exceeding {gap_threshold} seconds.")
            return list(zip(gaps.index, gaps.values))
        else:
            logger.info("No significant gaps detected in the data.")
            return []

    def save_dfdas_file(self, params=None):
        """
        Saves self.dfdas to an HDF5 file using self.filemeta.
        Automatically constructs the path if 'filename' is not in params.
        """
        if self.dfdas is None:
            logger.warning("No data in memory (self.dfdas) to save.")
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
            filename_out = save_to_DAS_file(
                self.dfdas, 
                filename, 
                "processed", 
                **self.filemeta
            )
            logger.info(f"Data successfully saved to: {filename_out[0]}")
        except Exception as e:
            logger.error(f"Error saving DFDAS file: {e}")


    def extract_channel_subset(self, offset_range_key):
        """
        Extracts a subset of channels from self.dfdas based on provided indices.
        Updates self.dfdas in place.
        """

        logger.info(f"Extracting channel subset for offset range key: {offset_range_key}")

        if offset_range_key not in self.offset_ranges:
            logger.warning(f"Offset range key '{offset_range_key}' not found in offset_ranges.")
            logger.warning(f"No changes made to dfdas.")
            logger.warning(f"Offset range key '{offset_range_key}' not found in offset_ranges.")
            return

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
            logger.warning("No data loaded in self.dfdas. Cannot extract channels.")
            return

        if self.dfdas_offsets is None:
            self.dfdas_offsets=self.dfdas.columns.to_numpy()

        original_offsets=set(self.dfdas_offsets)

        subset_dfdas = self.dfdas.columns.intersection(offsets)
        if subset_dfdas.empty:
            logger.warning("No matching channels found in self.dfdas for the provided offsets.")
            return

        missing_channels= set(offsets) - set(subset_dfdas)
        if missing_channels:
            # create a map of missing channels to their corresponding offset ranges for better clarity
            plot_channel_map = self.plot_channel_map(offsets, original_offsets)

            logger.warning(f"Warning: {len(missing_channels)} requested channels were not found in self.dfdas, and were skipped.")
            logger.info(f"Successfully extracted: {len(subset_dfdas)} channels.")
            logger.warning(f"The following requested channels were not found in self.dfdas: {sorted(list(missing_channels))}")
            logger.info(f"Successfully extracted: {len(subset_dfdas)} channels.")

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
            logger.warning(f"WARNING: 'channel_ranges' not specified for event '{event}'. Skipping channel calculation.")
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
            logger.info(f'Adjusting channels: {self.channels[0]} -> {self.channels[0] + adjustment}')
            self.channels = self.channels + adjustment
            
        return self.channels, self.min_channel, self.max_channel, self.channel_step

    def dfdas_to_csv(self, params=None):
        """
        Saves the decimated self.dfdas directly to CSV.
        Uses 'csv_out_path' and appends the timelapse suffix.
        """
        if self.dfdas is None:
            logger.warning("No data to save.")
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
        logger.info(f"»»» Saving CSV format to: {full_path}")
        # index=True includes the timestamps from the 't' axis created in decimation
        self.dfdas.to_csv(full_path, index=True, header=True)


    def dfdas_to_mseed(
        self, event: str, write: bool | str = False, format: str | None = None,
        dType: str | None = None, decimate: bool = False, filename: str | None = None,
        stationXML: bool = False, int_resolution: str = 'picostrain',
    ) -> Stream | None:
        """
        Converts the stored self.dfdas (Pandas) into an ObsPy Stream object.
        Uses class attributes for metadata and configuration.
        """
        # Short aliases for cleaner code
        ep, ap, pp, th = self.exec_params, self.acquisition_params, self.path_params, self.trace_header[self.network]

        if self.dfdas is None:
            logger.warning("No data loaded in self.dfdas. Run load_dfdas_file() first.")
            return None

        # reset decimation history for mseed output
        self._reset_decimation_history(domain='stream')

        use_int = dType == "int"
        if use_int:
            if int_resolution not in INT32_SCALES:
                raise ValueError(f"int_resolution must be one of {list(INT32_SCALES)}")
            scale = INT32_SCALES[int_resolution]
            logger.info(f"Scaling data by {scale:.0e} ({int_resolution}) for int32 MiniSEED.")
            # Stop here, before any trace is built, if the data won't fit.
            check_int32_scaling(self.dfdas.to_numpy(), scale)
            self.stationxml_params["int32_resolution"] = int_resolution
        else:
            # Float output: no counts mapping, so the key must not exist.
            # pop() with a default also clears a label left by an earlier int32 export.
            self.stationxml_params.pop("int32_resolution", None)

        if dType == 'int':
            logger.info(f'Scaling data by a factor of {scale:.0e} to match int32 range for MSEED output.\nNew data type is {dType}.')
        logger.debug(f'params: event:{event} write:{write} format:{format} dType:{dType} decimate:{decimate} filename:{filename}')
        if self.dfdas is None:
            logger.warning("No data loaded in self.dfdas. Run load_dfdas_file() first.")
            return None
        else:
            logger.info("Converting to ObsPy Stream")

        self.mseed_stream=Stream()
        params = self.get_current_params()   # built once, not once per channel
        for channel in self.dfdas.columns[:]:
            tr = Trace()
            data = self.dfdas[channel].to_numpy().astype(float)  # Convert to dType. default is float.
            tr.data = np.rint(data * scale).astype(np.int32) if use_int else data
            tr.stats.starttime = UTCDateTime(self.dfdas.meta["time"])
            tr.stats.sampling_rate = 1/self.dfdas.meta["dt"]

            station_code = str(int(channel*ap.get("dx"))).zfill(5)
            if ep.get("encode_station_name"):
                tr.stats['decoded_station'] = station_code
                tr.stats.station = encode_station_name(th.get("station_prefix"),station_code,params)
            else:
                tr.stats.station = station_code
            tr.stats.network = th['network']
            tr.stats.location = th['location']
            tr.stats.mseed = th['mseed']
            tr.stats.channel = th['channel']
            tr.stats.distance = int(channel*self.dfdas.meta["dx"])
            lat, long, ele = get_coordinates(channel, params, roi_channel=False)
            
            if ep.get("degrade_coordinates"):
                lat,long=degrade_coordinates(lat,long,precision=ep.get("precision"),error_metres=ep.get("error_metres"))
            
            tr.stats.coordinates = {'latitude': lat, 'longitude': long, 'elevation': ele}
            tr.stats.paz= AttribDict({'sensitivity': ap["sensitivity"]})
            # Append the trace to the Stream
            self.mseed_stream.append(tr)

        if decimate:
            self.mseed_stream=self.decimate_st()

        if stationXML:
            logger.debug("Creating StationXML file for the current stream.")
            station_xml_creation(self.mseed_stream,self.get_current_params())   

        if write in ['stream', True]:
            mseed_filename = filename if filename else f'{event}_{self.mseed_stream[0].stats.sampling_rate}Hz.mseed'
            full_path = os.path.join(pp.get("MSEED_path"), mseed_filename)
            self.mseed_stream.write(full_path, format='MSEED')
            logger.info(f'Full stream written to {full_path}')

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
                logger.info(f'wrote file {full_path} julian_date={ep.get("julian_date")}')

        elif not write:
            logger.info("Stream created in memory. No files written.")

        return self.mseed_stream
    
    def cross_correlate(self, reference=None, max_lag_seconds=10):
        """
        Cross-correlates all traces in self.dfdas and stores results in self.cross_correlation_results.
        """
        if self.dfdas is None:
            logger.warning("No data is loaded. Run load_dfdas_file() first.")
            return None
        else:
            logger.info("Performing cross-correlation on self.dfdas.")
        
        if reference is None:
            num_channels = len(self.dfdas.columns)
            reference_idx= (num_channels - 1) // 2
        else:
            reference= int(reference/self.roi_step) # Convert from channel number to index based on roi_step
            if reference >= len(self.dfdas.columns) or reference < -len(self.dfdas.columns):  # the negative index allows counting from the end, but it should still be within bounds
                logger.warning(f"Reference channel index {reference} is out of bounds. Defaulting to middle channel.")
                num_channels = len(self.dfdas.columns)
                reference_idx = (num_channels - 1) // 2
            else:
                reference_idx = reference

        reference_channel = self.dfdas.columns[reference_idx]
        channel_distance = int(reference_channel * self.dfdas.meta['dx'])
        logger.info(f"Using channel no. {reference_idx} as the reference channel for cross-correlation.")
        logger.info(f"Distance for reference channel: {channel_distance} metres.")

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
        logger.info(f"Cross-correlation completed with reference channel {reference_channel}.")
        return self.cross_correlation_results


    
    def decimate_st(self):
        """
        Decimates self.mseed_stream using the factors in freq_params.

        Each trace is filtered with the zero-phase Chebyshev II low-pass (filtfilt)
        and then decimated without a further filter. The effective zero-phase FIR
        of that filter is recorded on the stream (stream.decimation_steps), so the
        history lives exactly as long as the stream does.
        """
        fp = self.freq_params

        st = self.mseed_stream
        if st is None or len(st) == 0:
            logger.warning("No ObsPy Stream found in memory. Run dfdas_to_mseed() first.")
            return None

        for decimation in fp.get("decimation_factor"):
            fs = st[0].stats.sampling_rate
            logger.info(f'Decimation factor: {decimation}')
            new_nyquist = fs / 2.0 / decimation

            # Same coefficients zerophase_chebyshev_lowpass_filter will use, so
            # the recorded FIR is derived from exactly what is applied.
            b, a = get_cheby2_coeffs(fs * 0.5, new_nyquist)
            fir = zerophase_effective_fir(b, a)

            for tr in st:
                zerophase_chebyshev_lowpass_filter(tr, new_nyquist)
            st.decimate(factor=decimation, no_filter=True)

            # New list each time: never mutate a list another object may hold
            st.decimation_steps = list(getattr(st, "decimation_steps", [])) + [
                self._make_step(
                    decimation, fs, fir=fir,
                    label="Zero-phase Chebyshev II (filtfilt) low-pass")]

            logger.debug("%s", st[0].stats.sampling_rate)
        return st

    def reduce_velocity_st(self, vred=5.0, offset=0):
        """Return a copy of the MiniSEED stream with a linear moveout applied.

        Adds ``offset`` meters to each trace's distance, then shifts its
        start time earlier by that distance divided by ``vred`` (in km/s).
        The result is also stored in ``self.mseed_reduced_velocity``.

        Args:
            vred: Reduction velocity in kilometers per second.
            offset: Distance offset in meters applied to every trace.

        Returns:
            The reduced-velocity copy of ``self.mseed_stream``.
        """
        logger.info(f'Reducing velocity to {vred} km/s with an offset of {offset} meters.')

        st=self.mseed_stream.copy()
        
        if offset != 0:
            logger.info(f'An offset of {offset} meters with be added to every distance.')

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

        fp = self.freq_params

        if self.dfdas is None:
            logger.warning("No data loaded.")
            return False
        
        # If factors are passed via process_chain, use them; 
        # otherwise, fall back to global_params
        if factors is None:
            factors = fp.get("decimation_factor", [])
        else:
            fp["decimation_factor"] = factors  # Update global_params if factors are passed directly
            logger.info(f"Default decimation factors updated to: {factors}")
        
        for decimation in factors:
             # 1. Current dimensions come from the data itself, not from meta
             nt_in, nx = self.dfdas.shape
             fs_in = 1 / self.dfdas.meta['dt']

             # 2. Filter once, apply it, and record that same array
             fir = design_decimation_fir(decimation)
             sig_decimated = sps.resample_poly(
                 self.dfdas.to_numpy(), up=1, down=decimation, axis=0,
                 window=fir, padtype='edge'
             )

             # 3. New timing from the new dimensions.
             # resample_poly returns ceil(nt_in / decimation) samples, so the
             # sample interval is exactly dt * decimation.
             nt_out = sig_decimated.shape[0]
             dt_out = self.dfdas.meta['dt'] * decimation

             # Start from the data's own first timestamp: meta['time'] goes stale
             # after trim_DAS_data (truncate copies meta unchanged).
             if isinstance(self.dfdas.index, pd.DatetimeIndex):
                 tstart = self.dfdas.index[0].to_pydatetime(warn=False)
             else:
                 tstart = self.dfdas.meta['time']
             t = create_time_axis(tstart, nt_out, dt_out)

             # 4. Rebuild as a DASDataFrame. The constructor recomputes
             # dimensionSizes / dimensionRanges from the new shape.
             meta_out = self.dfdas.meta.copy()
             meta_out.update(dt=dt_out, time=tstart)

             self.dfdas = DASDataFrame(
                 sig_decimated,
                 index=t,
                 columns=self.dfdas.columns,
                 meta=meta_out
             )

             # Consistency check: metadata must match the new dimensions
             if tuple(self.dfdas.meta["dimensionSizes"]) != self.dfdas.shape:
                 logger.warning(
                     f"dimensionSizes {tuple(self.dfdas.meta['dimensionSizes'])} "
                     f"!= data shape {self.dfdas.shape} after decimation.")

             # 5. Record exactly what was applied (kept until new data is loaded)
             self.dfdas_decimation_steps.append(self._make_step(
                 decimation, fs_in, fir=fir,
                 label="Linear-phase FIR (Kaiser) low-pass"))

             logger.info(f"Decimated {nt_in} -> {nt_out} samples x {nx} channels; "
                   f"new sampling rate: {1/dt_out} Hz")

             # CRITICAL: Clear the cached ObsPy stream because the underlying
             # data (sampling rate/length) has changed. Its recorded steps
             # (stored on the stream) go with it.
             self.mseed_stream = None

        return self.dfdas
            

    def calculate_hilbert_envelope(self, params=None):
        """
        Calculates and saves the signal envelope plot.
        Supports overrides for cmap, vmin, and vmax via params.
        """
        if self.dfdas is None:
            logger.warning("No data in self.dfdas to plot.")
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

        logger.info(f"»»» Running envelope (Cmap: {cmap}, Vmin: {vmin}, Vmax: {vmax})")
        logger.info(f"»»» Running envelope (Cmap: {cmap}, Vmin: {vmin}, Vmax: {vmax})")
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


    def check_data_access(self, experiment_path, datetime_start):
        """
        Verifies that the experiment path and the date-specific 
        subfolder exist before attempting to read files.
        """
        logger.info(f"Checking data access in {experiment_path} for date {datetime_start.strftime('%Y-%m-%d')}")
        # 1. Check base path
        if not os.path.exists(experiment_path):
            logger.error(f"ERROR: Base path does not exist: {experiment_path}")
            return False

        # 2. Check date-specific subfolder (SimpleDAS standard: path/YYYYMMDD)
        date_str = datetime_start.strftime('%Y%m%d')
        date_path = os.path.join(experiment_path, date_str)
        
        if not os.path.exists(date_path):
            logger.error(f"ERROR: Date folder missing: {date_path}")
            logger.info("Check if the drive is mounted or if the event date is correct.")
            return False

        # 3. Check for data type subfolder (e.g., 'dphi')
        dphi_path = os.path.join(date_path, 'dphi')
        if not os.path.exists(dphi_path):
             logger.warning(f"WARNING: 'dphi' folder not found in {date_path}. SimpleDAS might fail.")
             return False
        
        return True
    
    # ============================================================
    # PLOT METHODS
    # ============================================================


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
            logger.error("ObsPy FDSN client not available. Install ObsPy to enable downloads.")
            return None

        # Default EIDA node (can be changed later)
        if eida_url is None:
            eida_url = 'https://webdc.eu'  # hardcoded EIDA-compatible node

        # Ask for station code if not provided
        if station_code is None:
            station_code = input("Enter station code (e.g. ABCD): ").strip()
            if not station_code:
                logger.error("No station code provided. Aborting download.")
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
                logger.error(f"Invalid start time: {e}")
                return None
        else:
            starttime = UTCDateTime(starttime)

        if endtime is None:
            e = input("End time (YYYY-MM-DDTHH:MM:SS): ").strip()
            try:
                endtime = UTCDateTime(e)
            except Exception as ex:
                logger.error(f"Invalid end time: {ex}")
                return None
        else:
            endtime = UTCDateTime(endtime)

        # Instantiate client and fetch waveforms
        try:
            client = Client(base_url=eida_url)
            logger.info(f"Connecting to EIDA node: {eida_url}")
            st = client.get_waveforms(network, station_code, location, channel, starttime, endtime)
            logger.info(f"Downloaded {len(st)} trace(s) for {station_code} from {starttime} to {endtime}.")
        except Exception as err:
            logger.error(f"Error downloading data from EIDA: {err}")
            return None

        # Store stream on instance
        if save_stream:
            if not hasattr(self, 'station_streams') or self.station_streams is None:
                self.station_streams = {}
            self.station_streams[station_code] = st

        return st

    def fk_filter(self, vmin=1400, vmax=None, sigma=(10, 10), taper=None,
                  dtype=np.float32):
        """
        Apparent-velocity (f-k) pass filter on self.dfdas, in place.

        Keeps energy with vmin <= |f/k| <= vmax (m/s); vmax=None means no upper
        limit. sigma is the Gaussian smoothing of the mask edge in
        (frequency bins, wavenumber bins).

        taper: optional Tukey window applied to the data before the FFT to reduce
            edge leakage. None disables it (default). A float alpha applies the
            same fraction to both axes; a tuple (alpha_time, alpha_space) sets
            them separately (None or 0 disables an axis). alpha=0 is
            rectangular, 1 is Hann; 0.05-0.2 is typical. The output is tapered too.
        dtype: working precision. float32 (default) halves memory and keeps the
            FFT in single precision (scipy.fft preserves it); pass np.float64
            for double precision.

        Memory: one real-FFT spectrum (about half the size of a full fft2) plus
        one mask; the mask comes from broadcasting f[:, None] against k[None, :],
        so no meshgrid arrays are built.
        """
        if self.dfdas is None:
            logger.warning("No data loaded; nothing to filter.")
            return None

        logger.info(f"Applying FK filter with vmin={vmin}, vmax={vmax}, taper={taper}.")

        dt = self.dfdas.meta['dt']
        dx = self.dfdas.meta['dx']

        # Copy only if the taper will modify the array in place
        data = self.dfdas.to_numpy().astype(dtype, copy=taper is not None)
        nt, nx = data.shape                 # actual shape, not stale meta

        if taper is not None:
            a_t, a_x = (taper, taper) if np.isscalar(taper) else taper
            for a in (a_t, a_x):
                if a is not None and not 0 <= a <= 1:
                    raise ValueError(f"Tukey alpha must be in [0, 1], got {a}.")
            if a_t:
                data *= sps.windows.tukey(nt, alpha=a_t).astype(dtype)[:, None]
            if a_x:
                data *= sps.windows.tukey(nx, alpha=a_x).astype(dtype)[None, :]

        # Real FFT along time (axis 0), full FFT along space (axis 1).
        # axes=(1, 0): the real transform runs over the last listed axis.
        spec = sfft.rfft2(data, axes=(1, 0), workers=-1)       # (nt//2+1, nx)
        del data

        f = sfft.rfftfreq(nt, d=dt).astype(dtype)[:, None]      # (nt//2+1, 1)
        k = np.abs(sfft.fftfreq(nx, d=dx)).astype(dtype)[None, :]  # (1, nx)

        # Velocity test without dividing by k: |f/k| >= vmin  <=>  f >= vmin*k
        mask = f >= vmin * k
        if vmax is not None:
            mask &= f <= vmax * k

        # Smooth the mask edges to reduce ringing. The mask depends on |f|, so a
        # blur over the full f-k plane would see values mirrored about f=0:
        # "mirror" on the frequency axis reproduces that exactly (even nt; for
        # odd nt the top rows differ by <3e-3). The wavenumber axis is periodic.
        mask = ndimage.gaussian_filter(mask.astype(dtype), sigma=sigma,
                                       mode=("mirror", "wrap"))
        spec *= mask
        del mask

        out = sfft.irfft2(spec, s=(nx, nt), axes=(1, 0), workers=-1)   # (nt, nx), real

        self.dfdas = DASDataFrame(
            out,
            index=self.dfdas.index,
            columns=self.dfdas.columns,
            meta=self.dfdas.meta.copy(),
        )
        self.mseed_stream = None   # cached stream no longer matches the data
        return self.dfdas


    def _resolve_whale_stream(
        self, event: str | Stream | Trace | None, stream: Stream | Trace | None
    ) -> Stream | Trace:
        """Pick the data to analyse, without touching self.mseed_stream when the
        caller supplies one.

        Priority: an explicit ``stream`` (or a Stream/Trace passed as ``event``),
        then the already-loaded self.mseed_stream, then a new stream built from
        ``event``. Validation of the stream itself is done by the detector.
        """
        if stream is None and isinstance(event, (Stream, Trace)):
            stream, event = event, None
        if stream is not None:
            return stream

        if self.mseed_stream is None or len(self.mseed_stream) == 0:
            if event is None:
                raise ValueError(
                    "No data to analyse: pass a Stream (stream=...) or an event key.")
            logger.info("Stream is empty; building it from event %s.", event)
            self.mseed_stream = self.dfdas_to_mseed(event, write=False)
            if self.mseed_stream is None or len(self.mseed_stream) == 0:
                raise ValueError(f"Could not build a stream for event {event!r}.")
        return self.mseed_stream


    def baleen_whale_detector(
        self, event: str | None = None, species: str | None = None,
        template: str | None = None, template_fs: float = 100,
        detection_mode: str = 'single', channel_idx: int | None = None,
        band_detect: tuple[float, float] = (15, 25), env_win: float = 0.2,
        thresh: float = 4.5, min_dist: float = 8.0, search_win_sec: float = 1.0,
        corr_thresh: float = 0.15, spec_source: str = "raw",
        spec_fmin: float = 10, spec_fmax: float = 40, spec_nperseg: int = 1042,
        spec_noverlap: int = 900, spec_nfft: int = 2048, spec_dr: float = 15,
        cmap: str = "jet", stream: Stream | Trace | None = None,
        report: bool = True, plot: bool = True,
        save_path: str | os.PathLike | None = "Whale_detector.png", show: bool = True,
    ) -> list[dict[str, Any]]:
        """
        Detect template-matched baleen whale calls in DAS data.

        When ``template`` is omitted, the bundled whale templates are initialized
        in the current directory and ``species`` selects the template.
        The detector selects one channel or processes all channels, finds candidate
        signal peaks from a smoothed amplitude envelope, and checks candidates by
        cross-correlating them with a whale-call template. It prints detection
        summaries and statistics, then displays and saves a plot for one processed
        channel as ``Whale_detector.png``.

        The data source is, in order of priority: ``stream`` (an ObsPy
        ``Stream``, or a single ``Trace``), the processor's already-loaded
        MiniSEED stream, or a stream built from ``event``. A supplied
        ``stream`` is analysed as-is and is neither modified nor stored on the
        processor. Traces may have different sampling rates; each is processed
        at its own rate.

        Args:
            event: Event key used to create the MiniSEED stream when no
                   ``stream`` is given and none is loaded. A Stream/Trace
                   passed here is treated as ``stream``.
            species: Whale species used to select a default template when
                     ``template`` is not provided. Supported values are
                     ``"blue"``, ``"fin"``, ``"humpback"``, ``"sei"``,
                     and ``"minke"``.
            template: Explicit template filename. Takes precedence over
                      ``species``.
            detection_mode: Channel selection mode: ``"single"`` uses
                            ``channel_idx``, ``"auto"`` chooses the
                            channel with the strongest band-passed
                            signal, and ``"all"`` processes every
                            channel.
            channel_idx: Zero-based stream index required when
                         ``detection_mode="single"``.
            band_detect: Low and high cutoff frequencies, in Hz, for the
                         detection band.
            env_win: Moving-average envelope window duration, in seconds.
            thresh: Candidate amplitude threshold as a multiple of the median
                    envelope amplitude.
            min_dist: Minimum separation between candidate peaks, in seconds.
            template_fs: Sampling rate, in Hz, of the template.
            search_win_sec: Cross-correlation search window, in seconds.
            corr_thresh: Minimum cross-correlation value for accepting a
                         candidate.
            spec_source: Signal used to calculate the spectrogram (raw or
                         detection-filtered data).
            spec_fmin: Lower spectrogram frequency limit, in Hz.
            spec_fmax: Upper spectrogram frequency limit, in Hz.
            spec_nperseg: Spectrogram window length in samples.
            spec_noverlap: Number of overlapping samples between spectrogram
                           windows.
            spec_nfft: FFT length used for the spectrogram.
            spec_dr: Spectrogram dynamic range, in dB.
            cmap: Matplotlib colormap name for the spectrogram.
            stream: ObsPy ``Stream`` (or ``Trace``) to analyse directly.
                    Must not contain masked (gappy) data.
            report: Print the detection summary to stdout.
            plot: Draw the detection figure for one processed channel. The
                  figure is also kept as ``self.whale_figure``.
            save_path: File the figure is saved to; ``None`` to skip saving.
            show: Call ``plt.show()``. The figure is closed afterwards either
                  way, so batch runs do not accumulate open figures; the
                  ``self.whale_figure`` object can still be saved or reshown.

        Returns:
            A list of dictionaries, ``all_results``, containing one record for
            each accepted detection across the processed channel(s). It is empty
            when no detections pass the correlation threshold. Each record has:

            - ``sample`` (int): Detection start index in the channel trace.
            - ``time_s`` (float): Detection start time relative to the trace start,
                                  in seconds.
            - ``datetime`` (datetime.datetime): Absolute detection start time.
            - ``corr`` (float): Maximum template cross-correlation value.
            - ``shift_s`` (float): Correlation shift from the approximate start,
                                   in seconds.
            - ``amp_max`` (float): Maximum absolute amplitude in the detection
                                   window.
            - ``signal_rms`` (float): RMS amplitude in the detection window.
            - ``noise_rms`` (float): RMS amplitude in the preceding noise window.
            - ``snr`` (float): Ratio of signal RMS to noise RMS; NaN when the
                               noise RMS is zero.
            - ``channel_idx`` (int): Zero-based index of the channel.
            - ``channel_label`` (str): Channel station label, or a generated
              ``channel_<index>`` label when no station label is available.

            ``single`` and ``auto`` modes return detections for one channel;
            ``all`` mode combines detections from all channels. Printing, plotting
            and saving are optional (``report``, ``plot``, ``save_path``,
            ``show``); with all of them off the method is pure computation.

        Raises:
            ValueError: If neither ``species`` nor ``template`` is provided, if
                the species or detection mode is unsupported, if single mode
                is selected without a channel index, if no stream, loaded data
                or event is available, or if the stream is empty or gappy.
            TypeError: If ``stream`` is not a Stream or Trace.
            FileNotFoundError: If the selected template file does not exist.
            IndexError: If ``channel_idx`` is outside the stream.
        """
        st = self._resolve_whale_stream(event, stream)
        result = detect_baleen_whales(
            st, species=species, template=template, template_fs=template_fs,
            detection_mode=detection_mode, channel_idx=channel_idx,
            params=DetectionParams(
                band_detect=band_detect, env_win=env_win, thresh=thresh,
                min_dist=min_dist, search_win_sec=search_win_sec, corr_thresh=corr_thresh),
            spectrogram=SpectrogramParams(
                source=spec_source, fmin=spec_fmin, fmax=spec_fmax, nperseg=spec_nperseg,
                noverlap=spec_noverlap, nfft=spec_nfft, dr=spec_dr, cmap=cmap),
            report=report, plot=plot, save_path=save_path, show=show)
        self.whale_figure = result.figure
        return result.detections
