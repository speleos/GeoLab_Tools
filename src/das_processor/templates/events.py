# ECHO_EQ-sources.py
# contains earthquake source parameters for teleseismic events M>5.0 , and regional events M>2.0


channel_ranges={
    "broadband_test_zoom": {"channel_min": [7000] ,
                "channel_max": [8500] ,    #[11293] ,
                "channel_step": [1]},#10

    "seismic_broadband_test": {"channel_min": [5000] ,
                "channel_max": [8500] ,    #[11293] ,
                "channel_step": [1]},#10

    "seismic": {"channel_min": [450] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [1]},#10

    "cross_correlation": {"channel_min": [3000] ,
                "channel_max": [11000] ,    #[11293] ,
                "channel_step": [5]},#10
    
    "seismic5": {"channel_min": [1000] ,
                "channel_max": [11283] ,    #[11293] ,
                "channel_step": [5]},


    "seismic_spectrum": {"channel_min": [450] ,
                        "channel_max": [4000] ,    #[11293] ,
                        "channel_step": [1]},#10
    
    "lowpass": {"channel_min": [3000] ,
                "channel_max": [10000] ,    #[11293] ,
                "channel_step": [100]},#10
    
    "internal_tides": {"channel_min": [500] ,
                "channel_max": [11000] ,    #[11293] ,
                "channel_step": [100]},#10
    
    "seismic_zoom": {"channel_min": [6500] ,
                "channel_max": [7000] ,
                "channel_step": [10]},

    "seedlink": {"channel_min": [392] ,
                "channel_max": [11283] ,
                "channel_step": [392]},

    "matthew_internal_waves": {"channel_min": [1125] ,
                "channel_max": [3083] ,
                "channel_step": [1]},

    "MUSIC": {"channel_min": [2000,3000,4000] ,
              "channel_max": [2100,3100,4100] ,
              "channel_step": [10,10,5]},
    "seismic_BEAM": {"channel_min": [4000] ,
             "channel_max": [4500] ,
             "channel_step": [100]},
    "nsd": {"channel_min": [0],
             "channel_max": [450],
             "channel_step": [1]},
    "seismic_BEAM_multiple": {"channel_min": [5250,1000,2000,3000,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,10500,11000,5640],
                      "channel_max": [5750,1500,2500,3500,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,10500,11000,11200,5680],
                      "channel_step": [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,1]},
    "coherence": {"channel_min": [450,650,850,1050,1250,1450,1650,1850,2050,2250,2450,2650,2850,3050,3250,3450,3650,3850,4050,4250,4450,
                                  4650,4850,5050,5250,5450,5650,5850,6050,6250,6450,6650,6850,7050,7250,7450,7650,7850,8050,8250,8450,
                                  8650,8850,9050,9250,9450,9650,9850,10050,10250,10450,10650,10850,11050],
                  "channel_max": [650,850,1050,1250,1450,1650,1850,2050,2250,2450,2650,2850,3050,3250,3450,3650,3850,4050,4250,4450,4650,
                                  4850,5050,5250,5450,5650,5850,6050,6250,6450,6650,6850,7050,7250,7450,7650,7850,8050,8250,8450,8650,
                                  8850,9050,9250,9450,9650,9850,10050,10250,10450,10650,10850,11050,11250],
                  "channel_step": [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]},
    "coastal": {"channel_min": [450] ,
                "channel_max": [2500] ,
                "channel_step": [1]},
    "madalia_world": {"channel_min": [300] ,
                "channel_max": [2500] ,
                "channel_step": [1]},
    "land": {"channel_min": [0] ,
              "channel_max": [800] ,
              "channel_step": [1]},
    "tecva": {"channel_min": [500] ,
              "channel_max": [11293] ,
              "channel_step": [500]},
    "spectrogram": {"channel_min": [200,1200,2200,3200,4200,5200,6200,7200,8200,9200,10200,11200] ,
                "channel_max":     [305,1205,2205,3205,4205,5205,6205,7205,8205,9205,10205,11205] ,    #[11293] ,
                "channel_step":    [10,  1,   1,1,1,1,1,1,1,1,1,1]},#10
    "spectrogram_currents": {"channel_min": [3000] ,
                "channel_max": [3005] ,    #[11293] ,
                "channel_step": [1]},#10
    "whales": {"channel_min": [8204] ,
                "channel_max": [8204] ,    #[11293] ,
                "channel_step": [1]},#10
    "blue_whale": {"channel_min": [3400] ,
                "channel_max": [3600] ,    #[11293] ,
                "channel_step": [20]},#10
    "ship_scrub": {"channel_min": [700] ,
                "channel_max": [1000] ,    #[11293] ,
                "channel_step": [5]},#10
    "whale_envelopes": {"channel_min": [2000] ,
                "channel_max": [8000] ,    #[11293] ,
                "channel_step": [1]},#10
    "whale_envelopes_R1": {"channel_min": [2000] ,
                "channel_max": [4000] ,    #[11293] ,
                "channel_step": [1]},#10
    "whales_zoom": {"channel_min": [8190] ,
                "channel_max": [8215] ,    #[11293] ,
                "channel_step": [1]},#10
    "snr_calculation": {"channel_min": [700] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [1]},#10
    "mseed_EIDA": {"channel_min": [200] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [800]},#10
    "mseed_IRIS": {"channel_min": [200] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [20]},#10
    "mseed_IRIS_1": {"channel_min": [200] ,
                "channel_max": [4000] ,    #[11293] ,
                "channel_step": [20]},#10
    "mseed_IRIS_2": {"channel_min": [4020] ,
                "channel_max": [8000] ,    #[11293] ,
                "channel_step": [20]},#10
    "mseed_IRIS_3": {"channel_min": [8020] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [20]},#10
    "mseed_IRIS_ALL": {"channel_min": [200] ,
                "channel_max": [11293] ,    #[11293] ,
                "channel_step": [20]},#10
}

sources={ 
"cross_correlation_test":
	{"start_time":'2026-04-12 08:30:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'cross_correlation' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'cross_correlation',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "save_DAS_event": False,
	"mode": 'seismic' ,
    "end_time":'2026-04-13 00:00:00.0',
	"window_duration": 14400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [10,10]}},
		{"operation":'cross_correlate', "params":{"reference": 7000, "max_lag_seconds": 1}},
		]
    },
         
"GEOFON_MSEED_20240804-08_A":
	{"start_time":'2024-08-03 08:30:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1', 
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240803_A":
	{"start_time": "2024-08-03T08:30:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1', 
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-04 00:00:00.0',
	"window_duration": 55800,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240803_B":
	{"start_time": "2024-08-03T08:30:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2', 
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-04 00:00:00.0',
	"window_duration": 55800,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240803_C":
	{"start_time": "2024-08-03T08:30:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3', 
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-04 00:00:00.0',
	"window_duration": 55800,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240804-08_B":
	{"start_time": "2024-08-04T00:00:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240804-08_C":
	{"start_time": "2024-08-04T00:00:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-08-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },

"GEOFON_MSEED_20260311-14_A":
	{"start_time":'2026-03-11 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-15 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20260311-14_B":
	{"start_time":'2026-03-11 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-15 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20260311-14_C":
	{"start_time":'2026-03-11 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-15 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },

"TESTE_INT":
	{"start_time":'2026-04-01 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-01 00:00:10.0',
	"window_duration": 10,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },


"GEOFON_MSEED_20260409-15_A":
	{"start_time":'2026-04-09 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260409-15_B":
	{"start_time":'2026-04-09 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260409-15_C":
	{"start_time":'2026-04-09 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-16 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },


"GEOFON_MSEED_20260401-08_A":
	{"start_time":'2026-04-03 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-09 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260401-08_B":
	{"start_time":'2026-04-08 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-09 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260401-08_C":
	{"start_time":'2026-04-01 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi/', "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-09 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
    
"GEOFON_MSEED_20260328-31_A":
	{"start_time":'2026-03-28 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260328-31_B":
	{"start_time":'2026-03-28 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
"GEOFON_MSEED_20260328-31_C":
	{"start_time":'2026-03-28 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels', "dType": 'int'}},
		]
    },
    
"GEOFON_MSEED_20260315-24_A":
	{"start_time":'2026-03-15 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-25 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20260315-24_B":
	{"start_time":'2026-03-15 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-25 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20260315-24_C":
	{"start_time":'2026-03-15 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-03-25 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },

"GEOFON_MSEED_20241203_13_A":
	{"start_time":'2024-12-03 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/mnt/SUBMERSE/DATA_dphi/' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-12-14 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20241203_13_B":
	{"start_time":'2024-12-03 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/mnt/SUBMERSE/DATA_dphi/' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-12-14 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20241203_13_C":
	{"start_time":'2024-12-03 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/mnt/SUBMERSE/DATA_dphi/' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-12-14 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },

"GEOFON_MSEED_20240816_31_A":
	{"start_time":'2024-08-16 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-09-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240816_31_B":
	{"start_time":'2024-08-16 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_2',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-09-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
"GEOFON_MSEED_20240816_31_C":
	{"start_time":'2024-08-16 00:00:00.0', "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_3',
    "DAS_data_path": '/media/loureiro/SUBMERSE 4/DATA_dphi' , "mseed_output_path": '/mnt/DOMLUIZ/02_MSEED/daily_mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2024-09-01 00:00:00.0',
	"window_duration": 86400,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },
    
"TEST-S3":
	{"start_time": "2023-10-29 19:09:25" ,
	"end_time": "2023-10-29 19:09:45" ,
    "mode": 's3_download' , # 'channel_ranges': 'mseed_IRIS_1', #doesn't matter for s3 download, but needed for the process chain
    "DAS_data_path": '/tmp/S3-test/', 
    "s3_bucket": 'gc.3x2023', "endpoint_url": 'https://s3.gfz-potsdam.de/',"unsigned": True,} ,

"TEST-MSEED":
	{"start_time": "2023-10-29 19:09:25" ,
	"end_time": "2023-10-29 19:09:45" ,
    "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_1',
    "mode": 'seismic' ,
    "DAS_data_path": '/tmp/S3-test/', "mseed_output_path": '/tmp/S3-test/',
    'process_chain':[
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'stream'}},
		]
    } ,

"STATIONXML_20260413":
	{"start_time": "2026-04-13T00:00:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_ALL',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/tmp/mseed' ,
	"mode": 'time_lapse' ,
    "end_time":'2026-04-13 00:00:30.0',
	"window_duration": 30,
    'process_chain':[
		{"operation":'decimate_dfdas', "params":{"factors": [2]}},
        {"operation":'dfdas_to_mseed', "params":{"write": 'channels'}},
		]
    },

"cut_test":
	{"start_time": "2026-04-13T00:00:00" , "magnitude": 2.3   , "latitude": 31.51 , "longitude": -16.21 , "depth": 0.0 ,
	"type":'seismic' , "window_in": 0 ,  "local": True , "GPS": True, "channel_ranges": 'mseed_IRIS_ALL',
    "DAS_data_path": '/mnt/DOMLUIZ/DATA_dphi' , "mseed_output_path": '/tmp/mseed' ,
	"mode": 'seismic' ,
    "end_time":'2026-04-13 00:00:30.0',
	"window_duration": 30,
    'process_chain':[
		{"operation":'extract_channel_subset', "params":{"offset_range_key": "test"}},
		]
    },


}

sources_test = {
	"version": 2.0 ,
    
"2023OCT29_1909_TEST-S3":
	{"start_time": "2023-10-29 19:09:25" ,
	"end_time": "2023-10-29 19:09:45" ,
    "mode": 's3_download' ,
    "DAS_data_path": '/tmp/S3-test/', 
    "s3_bucket": 'gc.3x2023', "endpoint_url": 'https://s3.gfz-potsdam.de/',"unsigned": True,} ,
}

