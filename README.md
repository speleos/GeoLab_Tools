# GeoLab tools
Tools to automate the processing of OptoDAS datasets
Copyright (C) 2026 Afonso Loureiro / ARDITI

## Usage

There are two configuration files:

- *defaults.py*, where general parameters that modify function behaviour are defined.
- *sources.py*, where the events and processing steps are defined

A commented example entry for *sources.py*:

```
sources = {
    # Name of the event. It is also used for filenames of saved data and figures.
    "MSEED_20231028-30_A": {

        # The event time
        "start_time": '2023-10-28 00:00:00.0',

        # Only needed if running time lapses
        "end_time": '2023-10-31 00:00:00.0',

        # Used when processing earthquake data
        "magnitude": 2.3,

        # Coordinates are used to calculate theoretical travel times
        "latitude": 31.51,
        "longitude": -16.21,
        "depth": 0.0,

        # Processing type: If 'seismic' the processor assumes that start_time is event time
        # and will calculate the travel time to select the appropriate data
        # Other defined types are 'spectrogram', 'whales', 'spectrogram_currents'
        "type": 'seismic',

        # Available modes are 'time_lapse', which applies the processing chain as many times
        # as needed between 'start_time' and 'end_time', and 'seismic', which uses 'start_time'
        # and 'window_duration'
        "mode": 'time_lapse',

        # Defines how much time before the event should be loaded
        "window_in": 0,

        # How much data should be loaded
        "window_duration": 3600,

        # Set to True for local and regional seismic events
        "local": True,

        # Indicate if precise timing is available, and apply a clock correction if needed.
        # It will be checked upon data loading.
        "GPS": True,

        # Which channels to load, taken from a dictionary called channel_ranges
        "channel_ranges": 'mseed_A',

        # Where the raw data is to be found or saved.
        # This overrides any path definitions in shared_config.py
        "DAS_data_path": '/mnt/data/S3-test',

        # Where miniSEED files are to be saved. Also overrides previously-defined paths.
        "mseed_output_path": '/mnt/data/02_MSEED/hourly_mseed',

        # The operations to perform after loading the data, in the order they should be performed.
        "process_chain": [
            {
                # The available operations are the functions in processor.py
                "operation": 'decimate_dfdas',
                # The parameters are also the ones defined for each function.
                # Parameters will be matched to the function.
                "params": {"factors": [2]}
            },
            {
                # Multiple invocations of the same function are possible
                "operation": 'dfdas_to_mseed',
                # in this case, the first invocation produces a single miniSEED
                # file with as many traces as loaded DAS channels
                "params": {"write": 'stream'}
            },
            {
                "operation": 'dfdas_to_mseed',
                # in the second invocation, one file per channel will be saved
                # if 'write' is set to false, no files are written, and only the OBSpy Stream() is created
                "params": {"write": 'channels'}
            },
        ]
    },
}
```
---
## License & Attribution

This project is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)** - see the [LICENSE](LICENSE) file for details.

* **Copyright © 2026 Afonso Loureiro / ARDITI** - Tools and Framework
* **Portions Copyright © ASN** - Certain algorithms, modules, or foundational components are copyright **ASN** and used under the terms of the GNU General Public License (GPL-3.0-or-later).
---
## Funding
This development was supported by the European Commission under the Horizon Europe programme, through the project SUBMERSE (https://submerse.eu/) (Grant Agreement No. 101095055), funded within the HORIZON-INFRA-2022-TECH-01 call.
Additional support from the ECHO project (DOI: 10.54499/2024.13655.PEX), and FCT, I.P./MCTES through national funds (PIDDAC) (LA/P/0068/2020 - DOI: 10.54499/LA/P/0068/2020, UID/50019/2025 - DOI: 10.54499/UID/50019/2025, UID/PRR/50019/2025 - DOI: 10.54499/UID/PRR/50019/2025,  UID/PRR2/50019/2025 - DOI: 10.54499/UID/PRR2/50019/2025).
