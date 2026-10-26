
# default paths
path_params={
    "coords_path": '/mnt/D2A/DAS/ELLALink/',
    "csv_out_path": '/mnt/SUBMERSE/ECHO/CSV/',
    "csv_SNR_path": '/mnt/D2A/DAS/SNR/',
    "coords_DAS": 'roi_chan_lat_long_depth_dist_azim_dip_5m_17OCT2024.csv',
    "coords_DAS_lowres": 'ELLALink_public_coords.csv',
    "csv_RMS": 'channels_RMS.csv',
    "MSEED_path": '/mnt/D2A/DAS/MSEED/tmp',

    "SGY_path": '/mnt/D2A/DAS/SGY',
    "pickle_path": '/mnt/D2A/DAS/pickles',
    "DAS_data_path": '/mnt/INESCTEC/',
    
    "experiment_path": '/mnt/INESCTEC/',

    "events_path": '/mnt/SUBMERSE/Events/ECHO/',
    
    "mseed_land_stations": '/mnt/D2A/DAS/CEIDA/27OCT2023_all-stations.mseed',
    "xdas_path": '/mnt/SUBMERSE/XDAS',
    "hdf5_headers_example": '/mnt/D2A/DAS/OptoDAS_SN44/26OCT-02NOV/20231027/dphi/125959.hdf5',
    "JSON_template": '/mnt/D2A/DAS/DAS-RCN_example.json',
    "figures_path": '/mnt/D2A/DAS/figures/'
}

clock_offset=150

# default parameters for frequency analysis
freq_params={
    "decimation_factor": [8], #[5]
    "land_station_min_freq": 1.,
    "land_station_max_freq": 50.,
}

station_coordinates = {
    #"PMAR": {"latitude": 32.72, "longitude": -16.91, "channels": "SH*"},
    "PMAR": {"latitude": 32.72, "longitude": -16.91, "channels": "HH*"},
    "PMPST": {"latitude": 33.08, "longitude": -16.33, "channels": "HH*"},
   # "PMPST": {"latitude": 33.08, "longitude": -16.33, "channels": "HN*"},
    "PMPS": {"latitude": 33.06, "longitude": -16.33, "channels": "HN*"},
    "PMOZ": {"latitude": 32.82, "longitude": -17.20, "channels": "HN*"},
    "FUL": {"latitude": 32.65, "longitude": -16.89, "channels": "SH*"},
    "FIMU":{"latitude": 32.29856, "longitude": -16.62694},
}

# default acquisition parameters for OptoDAS
acquisition_params={
    "roi_step": 5, #defaults for OptoDAS
    "gauge_length": 10,
    "dx": 1,
    "overlap": 0.5,
    "gps_status": 6, #default is 6=no GPS. 0=GPS locked.
    "sensitivity": 1,
    "sensitivityUnits": "rad/(strain*m)",
    "original_sample_rate": 500,
}

trace_header={'5V':
              {
                'network': '5V',
                "location": '',
                "mseed":'D',
                "channel":'HSF',
                "net_description": "GeoLAB DAS cable. 57 km-long dark fibre starting at the Praia Formosa CLS, managed by ARDITI.",
                "net_starttime": [2023,10,26],
                "net_no_stations": 11293, # load this value later
                "sta_name_prefix": "GeoLAB fibre. DAS",
                "inventory_source": "maloureiro@fc.ul.pt",
                "station_prefix": "G", #"GLB", # the number of letters in the station name prefix will change the precision of the number code.
             },
             '3X':
             {
                "network": '3X',
                "location": '',
                "mseed":'D',
                "channel":'HSF',
                "net_description": "GeoLAB DAS cable. 57 km-long dark fibre starting at the Praia Formosa CLS, managed by ARDITI.",
                "net_starttime": [2023,10,26],
                "net_endtime": [2023,11,3],
                "net_no_stations": 11293, # load this value later
                "sta_name_prefix": "GeoLAB fibre. DAS",
                "inventory_source": "maloureiro@fc.ul.pt",
                "station_prefix": "G", #"GLB", # the number of letters in the station name prefix will change the precision of the number code.
             }
}

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

offset_ranges={
    "test": {"offset_min": [2300,10000] ,
             "offset_max": [8000,15000] ,  
                "offset_step": [5,5]},
}
