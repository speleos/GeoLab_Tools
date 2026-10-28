# das-processor

Tools for processing Distributed Acoustic Sensing (DAS) data recorded with OptoDAS interrogators: cut an event window out of continuous HDF5 files, decimate it, export it to MiniSEED / CSV / HDF5, build StationXML, and produce record sections, spectrograms and other plots.

Work is organised around **events**. An event is a dictionary in `sources.py` that says *what* to cut (time window, channels, location) and, optionally, *which processing steps* to run on it (the **process chain**). You can also skip the process chain and drive the processor yourself, step by step (see [Processing without a process chain](#processing-without-a-process-chain)).

## Contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [The configuration files](#the-configuration-files)
4. [Defining an event with processing steps](#defining-an-event-with-processing-steps)
5. [Running events](#running-events)
6. [Other event modes](#other-event-modes)
7. [Processing without a process chain](#processing-without-a-process-chain)
8. [Output files and logs](#output-files-and-logs)
9. [Troubleshooting](#troubleshooting)
10. [License and credits](#license-and-credits)

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/YOUR_USER/YOUR_REPO.git
```

Optional extras:

```bash
pip install "das-processor[s3] @ git+https://github.com/YOUR_USER/YOUR_REPO.git"        # S3 downloads (boto3)
pip install "das-processor[notebook] @ git+https://github.com/YOUR_USER/YOUR_REPO.git"  # progress bars in Jupyter
```

For development:

```bash
git clone https://github.com/YOUR_USER/YOUR_REPO.git
cd YOUR_REPO
pip install -e ".[dev]"
```

## Quick start

**1. Create the configuration files** in the folder where you will work:

```bash
das-processor init
```

or from Python:

```python
from das_processor import init_config
init_config()          # creates defaults.py and sources.py if they do not exist
```

Existing files are never overwritten (use `init_config(overwrite=True)` to reset them).

**2. Edit them.** In `defaults.py` set your data paths, channel ranges and station coordinates. In `sources.py` define the events you want to process.

**3. Run** (a complete script is shown in [Running events](#running-events)):

```python
from das_processor import *
from sources import sources
from defaults import *          # or import the names you need, see below

results = process_events(["MY_EVENT"], global_params)
DAS = results["MY_EVENT"]
DAS.dfdas                        # the loaded data
```

## The configuration files

`defaults.py` and `sources.py` are plain Python files that live in **your** working folder, not inside the package. Edit them freely.

### `defaults.py`

| Name | Contents |
|---|---|
| `path_params` | Folders for raw data (`experiment_path`, `DAS_data_path`), cut events (`events_path`), MiniSEED output (`MSEED_path`), figures (`figures_path`), coordinate tables and so on |
| `channel_ranges` | Named channel selections: lists of `channel_min`, `channel_max`, `channel_step` (in units of the ROI step). Events refer to them by name |
| `offset_ranges` | Named offset selections (used by `extract_channel_subset`) |
| `freq_params` | Default `decimation_factor` and land-station filter limits |
| `acquisition_params` | OptoDAS defaults: `roi_step`, `gauge_length`, `dx`, `overlap`, `gps_status`, sensitivity, `original_sample_rate`. Updated from the file headers while loading |
| `trace_header` | Network metadata per network code (`'5V'`, `'3X'`, ...) used for MiniSEED headers and StationXML |
| `station_coordinates` | Land/reference stations, with `latitude` and `longitude` (and `channels`) |
| `clock_offset` | Clock correction in seconds, applied to events where `GPS` is `False` |

### `sources.py`

Defines one dictionary called `sources`, with one entry per event. See the next section.

## Defining an event with processing steps

Add an entry to `sources` in `sources.py`. The key is the event name:

```python
sources = {
    "STATIONXML_20260413": {
        "start_time": "2026-04-13T00:00:00",
        "type": "seismic",
        "magnitude": 2.3, "latitude": 31.51, "longitude": -16.21, "depth": 0.0,
        "window_in": 0,
        "window_duration": 30,
        "GPS": True,
        "channel_ranges": "mseed_IRIS_ALL",
        "DAS_data_path": "/mnt/DOMLUIZ/DATA_dphi",
        "mseed_output_path": "/tmp/mseed",

        "process_chain": [
            {"operation": "decimate_dfdas", "params": {"factors": [2]}},
            {"operation": "dfdas_to_mseed", "params": {"write": "channels"}},
        ],
    },
}
```

### Event fields

| Field | Meaning |
|---|---|
| `start_time` | Reference time, ISO format. For `type: "seismic"` this is the earthquake origin time; for any other type it is the start of the window |
| `type` | `"seismic"` computes the P arrival with TauP at the reference station (needs `magnitude`, `latitude`, `longitude`, `depth`) and cuts the window around it. Any other value cuts `start_time - window_in` directly |
| `window_in` | Seconds to include **before** the P arrival (seismic) or before `start_time` (other types) |
| `window_duration` | Seconds to include after that point. Total length is `window_in + window_duration` |
| `GPS` | `True` if the interrogator clock was GPS-locked. If `False`, `clock_offset` is applied to the start time |
| `channel_ranges` | Name of an entry in `channel_ranges` in `defaults.py`. Validated at start-up: unknown names raise an error |
| `DAS_data_path` | Folder with the raw files (sub-folders named `YYYYMMDD`). Overrides `path_params["experiment_path"]` for this event only |
| `mseed_output_path` | Overrides `path_params["MSEED_path"]` for this event |
| `events_path` | Overrides `path_params["events_path"]` (where the cut event HDF5 is cached) |
| `save_DAS_event` | `True`/`False` overrides `exec_params["save_DFDAS"]` for this event |
| `mode` | Optional. See [Other event modes](#other-event-modes). Leave out for a normal event |
| `process_chain` | Optional list of processing steps (below) |

Path overrides are applied while the event runs and restored afterwards, so events with different data folders can be processed in one call.

### The process chain

`process_chain` is a list of steps run **in order** on the loaded data, each step acting on `DAS.dfdas` in place. Each step has the name of a `DASProcessor` method and its arguments:

```python
{"operation": "<method name>", "params": {"<argument>": <value>, ...}}
```

Rules:

- `operation` can be any method of `DASProcessor`. Unknown names are logged as a warning and skipped.
- `params` are passed as keyword arguments. Keys that the method does not accept are **silently dropped**, so check spelling.
- If the method has an `event` or `suffix` argument you do not need to supply it: the current event name and the window suffix are filled in automatically.
- Order matters. For example, decimate before converting to MiniSEED.

Commonly used operations:

| Operation | Main arguments | What it does |
|---|---|---|
| `decimate_dfdas` | `factors` (list, e.g. `[2]` or `[5, 2]`) | Anti-aliased decimation of `dfdas`. Without `factors`, uses `freq_params["decimation_factor"]` |
| `dfdas_to_mseed` | `write`, `format`, `dType`, `decimate`, `filename`, `stationXML` | Converts to an ObsPy `Stream` and optionally writes MiniSEED |
| `dfdas_to_csv` | `params` | Saves the decimated data to CSV |
| `save_dfdas_file` | `params` | Saves `dfdas` to HDF5 |
| `extract_channel_subset` | `offset_range_key` | Keeps the channels of a named entry in `offset_ranges` |
| `fk_filter` | `fs`, `vmin`, `vmax` | Frequency-wavenumber filter |
| `plot_seismic_section` / `plot_seismic_2d` / `plot_seismic_zoom` / `plot_seismic_enhanced` | `channel_step`, `start_offset_km`, `end_offset_km`, `max_time_s`, `scale_factor`, ... | Record sections (use travel-time phases, so they expect a `seismic` event) |
| `plot_spectrogram` | `freq_type` | Spectrogram plots |
| `calculate_hilbert_envelope` | `params` | Envelope calculation and plot |
| `cross_correlate`, `plot_cross_correlation_heatmap`, `plot_cross_correlation_lines` | `reference`, `max_lag_seconds`, ... | Cross-correlation of all channels and plots |
| `plot_record_sections` | `dataset` | Combined record sections (with land station) |
| `download_waveform_data` | `station_code`, `starttime`, `endtime`, `network`, `channel`, ... | Download land-station waveforms from an EIDA node |

Use `help(DASProcessor.<method>)` for the full list of arguments.

## Running events

Save as `run.py` in the folder that contains `defaults.py` and `sources.py`:

```python
from das_processor import *

from sources import sources
from defaults import (channel_ranges, offset_ranges, freq_params,
                      path_params, trace_header, station_coordinates,
                      clock_offset, acquisition_params)

events_to_process = ["STATIONXML_20260413"]

exec_params = {
    "events_to_process": events_to_process,
    "station": "FIMU",          # key of station_coordinates, used for epicentral distance
    "land_station": "PMAR",
    "network": "5V",            # key of trace_header
    "clock_offset": clock_offset,   # seconds, used when an event has "GPS": False

    "save_DFDAS": True,         # cache the cut window as <events_path>/<event>.hdf5
    "save_figures": True,
    "replace_figures": False,
    "stationXML": True,
    "encode_station_name": True,
    "julian_date": True,
    "degrade_coordinates": False,
    "verbose": True,
}

global_params = {
    "freq_params": freq_params,
    "exec_params": exec_params,
    "path_params": path_params,
    "trace_header": trace_header,
    "channel_ranges": channel_ranges,
    "offset_ranges": offset_ranges,
    "sources": sources,
    "station_coordinates": station_coordinates,
    "acquisition_params": acquisition_params,
}

results = process_events(events_to_process, global_params)

DAS = results["STATIONXML_20260413"]
print(DAS.dfdas.shape)        # (samples, channels)
```

What `process_events` does for every event in the list:

1. Looks up the event in `sources` (unknown names are reported and skipped).
2. Works out the time window (TauP travel times for seismic events, clock offset when `GPS` is `False`).
3. Reads the channels selected by `channel_ranges`, either from the cached event HDF5 in `events_path` (if present) or from the raw files in `DAS_data_path`.
4. Runs the event's `process_chain`, if it has one.
5. Writes a log file for the event.

`process_events` returns a dictionary `{event name: DASProcessor}`. The loaded data is in `results[event].dfdas`. Keep an eye on memory: all events stay in the dictionary until you delete them (`del results[event]`).

Run it with `python run.py`, or paste it into a notebook.

> `exec_params["station"]` must be a key of `station_coordinates`, and `exec_params["network"]` a key of `trace_header`.

## Other event modes

Set `"mode"` in the event to change what `process_events` does.

### `time_lapse`: many consecutive windows

Cuts the interval `start_time` to `end_time` into windows of `window_duration` seconds and runs the process chain on every window. Windows that already exist in `events_path` are skipped, and each window's outputs get a unique suffix such as `_20260413_0000`.

```python
"TIMELAPSE_20260413": {
    "mode": "time_lapse",
    "start_time": "2026-04-13T00:00:00",
    "end_time": "2026-04-13T01:00:00",
    "window_duration": 30,
    "window_in": 0,
    "type": "seismic", "GPS": True,
    "channel_ranges": "mseed_IRIS_ALL",
    "DAS_data_path": "/mnt/DOMLUIZ/DATA_dphi",
    "mseed_output_path": "/tmp/mseed",
    "process_chain": [
        {"operation": "decimate_dfdas", "params": {"factors": [2]}},
        {"operation": "dfdas_to_mseed", "params": {"write": "channels"}},
    ],
},
```

Data is released after each window, so `results[event].dfdas` is `None` for time-lapse events.

### `s3_download`: fetch raw files from S3

Needs the `s3` extra (`boto3`). Downloads the files covering `start_time` to `end_time` into `DAS_data_path`:

```python
"S3_20260413": {
    "mode": "s3_download",
    "s3_bucket": "my-bucket",
    "endpoint_url": "https://s3.example.org",
    "unsigned": True,                 # anonymous access
    "start_time": "2026-04-13T00:00:00",
    "end_time": "2026-04-13T01:00:00",
    "DAS_data_path": "/mnt/DOMLUIZ/DATA_dphi",
},
```

## Processing without a process chain

The process chain is a convenience. Everything it does is a normal method call, so you can load data and decide what to do next interactively. There are two levels.

### A. Load an event, then call the methods yourself

Define the event in `sources.py` **without** a `process_chain` (it still provides the time window, channels and paths), then work on the result:

```python
results = process_events(["STATIONXML_20260413"], global_params)
DAS = results["STATIONXML_20260413"]

DAS.dfdas                              # DASDataFrame: rows = time samples, columns = channels
DAS.dfdas.meta                         # metadata, e.g. meta["dt"] (sample interval) and meta["time"]

DAS.decimate_dfdas(factors=[2])        # same as the "decimate_dfdas" step
DAS.plot_seismic_section(channel_step=10, max_time_s=30)
DAS.dfdas_to_mseed("STATIONXML_20260413", write="channels")
```

You see the data between steps, repeat a step with other arguments, or stop after any step. The arguments are exactly the `params` you would put in a chain.

Without `process_events`, build the processor yourself:

```python
event = "STATIONXML_20260413"

DAS = DASProcessor(global_params)
DAS.event_info = sources[event]

with DAS.temporary_event_paths(event):      # applies the event's path overrides
    ok = DAS.load_event_data(event)         # True when data was found and loaded
    if ok:
        DAS.decimate_dfdas(factors=[2])
        DAS.dfdas_to_mseed(event, write="channels")
```

### B. Read the files directly, with no event definition

If you only want the data as a table, bypass `sources` completely and use the low-level readers:

```python
import datetime
import numpy as np
from das_processor import find_DAS_files, load_DAS_files

# 1. Locate the files overlapping a time window.
files, chIndex, samples = find_DAS_files(
    "/mnt/DOMLUIZ/DATA_dphi",                      # folder with YYYYMMDD sub-folders
    start=datetime.datetime(2026, 4, 13, 0, 0, 0),
    duration=datetime.timedelta(seconds=30),
    channels=np.arange(2000, 8000, 50),            # channel numbers as in the file header
)

# 2. Load them.
dfdas = load_DAS_files(files, chIndex, samples, integrate=False)

print(dfdas.shape)           # (samples, channels)
print(dfdas.meta["dt"])      # sample interval in seconds
```

`dfdas` is a pandas DataFrame, so you can use any scientific Python tool on it:

```python
import pandas as pd
from scipy import signal

fs = 1 / dfdas.meta["dt"]
sos = signal.butter(4, [1, 20], btype="bandpass", fs=fs, output="sos")

filtered = pd.DataFrame(
    signal.sosfiltfilt(sos, dfdas.to_numpy(), axis=0),
    index=dfdas.index, columns=dfdas.columns,
)
rms_per_channel = filtered.pow(2).mean().pow(0.5)
```

To reuse the package's own steps on this data, attach it to a processor. Steps that only need the data and its metadata, such as decimation, then work directly:

```python
DAS = DASProcessor(global_params)
DAS.dfdas = dfdas
DAS.decimate_dfdas(factors=[2])
```

Steps that need event information (station coordinates, travel times, MiniSEED/StationXML headers, plots with phases) also need `DAS.event_info`, `DAS.filemeta`, `DAS.datetime_start`, `DAS.duration` and `DAS.channels`. For those, use approach A.

### Other module-level tools

These do not need a processor or an event, and can be imported from `das_processor`:

| Function | Purpose |
|---|---|
| `get_data_indexes(path)`, `get_filemeta(path, metaDetail)` | Read channel/sample indexes and header metadata of one HDF5 file |
| `save_to_DAS_file(dfdas, filepath)` | Write a `DASDataFrame` to an HDF5 file |
| `get_travel_times(event, global_params, rec_lat, rec_lon)` | TauP P and S arrivals for an event defined in `sources` |
| `zerophase_chebyshev_lowpass_filter(trace, freqmax)` | Zero-phase Chebyshev II low-pass on an ObsPy trace |
| `degrade_coordinates(lat, lon, precision, error_metres)` | Blur coordinates for public data release |
| `baleen_whale_detector(st, global_params, ...)` | Fin-whale call detector on an ObsPy `Stream` |

## Output files and logs

| What | Where |
|---|---|
| Cut event window (cache) | `<events_path>/<event>.hdf5`, written when `save_DFDAS` is `True`. On later runs this file is loaded instead of the raw data. To force a reload from raw data, call `DAS.load_event_data(event, replace_DAS=True)` |
| MiniSEED | `path_params["MSEED_path"]`, or the event's `mseed_output_path` |
| Figures | `path_params["figures_path"]` |
| Log | A `.log` file in the current folder, named after the event and the start time |

## Troubleshooting

**`ModuleNotFoundError: No module named 'sources'`**: the configuration files are not in the current folder. Run `das-processor init` there, or start Python from the folder that contains them.

**`<event>: '<name>' does not exist in the defined channel_ranges`**: the event's `channel_ranges` value is not a key in `channel_ranges` in `defaults.py`.

**`KeyError` on `station_coordinates[station]`**: `exec_params["station"]` is not defined in `station_coordinates`.

**No files found / "Skipping event due to missing directories"**: the raw data folder must contain `YYYYMMDD` sub-folders. Check `DAS_data_path` in the event, or `experiment_path` in `path_params`, and that the dates cover the window.

**A step does nothing**: parameters that the method does not accept are dropped without warning. Run the step by hand (approach A) to see its signature, or use `help(DASProcessor.<method>)`.

## License and credits

Released under the GNU General Public License v3.0 (see `LICENSE`).

The DAS file reading code is derived from [simpleDAS](https://github.com/ASN-Norway/simpleDAS) by ASN (GPL-3.0). The package also builds on [ObsPy](https://docs.obspy.org/), NumPy, SciPy, pandas and Matplotlib.
