sources={ 
"cross_correlation": # performs cross-correlation between traces
	{"start_time":'2026-04-12 08:30:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'cross_correlation' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'cross_correlation',
    "DAS_data_path": './DATA2/DATA_dphi' , # overrides the default path for DAS data
	"save_DAS_event": False,
	"mode": 'seismic' ,
    "end_time":'2026-04-13 00:00:00.0',
	"window_duration": 14400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [10,10]}}, # decimating 100 times, in two steps
		{"operation":'cross_correlate', "params":{"reference": 7000, "max_lag_seconds": 1}},
		]
    },
         
"DECIMATE_MSEED": # loads a day of data, decimates it, saves individual channels as separate miniSEED files
	{"start_time":'2024-08-03 08:30:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_A', 
    "DAS_data_path": './DATA2/DATA_dphi' , # overrides the default path for DAS data
	"mseed_output_path": './02_MSEED/daily_mseed' , # overrides the default output path for mseed files
	"mode": 'time_lapse' ,
    "end_time":'2024-08-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
   
"S3 download": # downloads data from remote S3 storage
	{"start_time": "2023-10-29 19:09:25" ,
	"end_time": "2023-10-29 19:09:45" ,
    "mode": 's3_download' ,
    "DAS_data_path": './S3-DATA/', # overrides the default path for DAS data
    "s3_bucket": 'gc.3x2023', "endpoint_url": 'https://s3.gfz-potsdam.de/',"unsigned": True,} ,

"DECIMATE-MSEED-MSEED":  # loads a day of data, decimates it, saves individual channels as separate miniSEED files and a single file with all traces, will also create a stationXML file
	{"start_time": "2023-10-29 19:09:25" ,
	"end_time": "2023-10-29 19:09:45" ,
    "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_A',
    "mode": 'seismic' ,
    "DAS_data_path": '/tmp/S3-test/', "mseed_output_path": '/tmp/S3-test/',
    'process_chain':[
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'stream', "stationXML": True}},
		]
    } ,

"SUBSET": # if event data has already been saved with this name, it will load it and extract a sbuset of offsets
		  # if not, it will load the data according to the supplied channel ranges and then extract the subset
	{"start_time": "2026-04-13T00:00:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_ALL',
	"mode": 'seismic' ,
    "end_time":'2026-04-13 00:00:30.0',
	"window_duration": 30,
    'process_chain':[
		{"operation":'extract_channel_subset', "params":{"offset_range_key": "test"}},
		]
    },
}

