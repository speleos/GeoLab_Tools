
# default paths
path_params={
    "coords_path": './COORDS/',
    "coords_DAS": 'roi_chan_lat_long_depth_dist.csv',
    "coords_DAS_lowres": 'low_res_coords.csv',

    "csv_out_path": './CSV/',
    "csv_SNR_path": './SNR/',
    "csv_RMS": 'channels_RMS.csv',
    
    "MSEED_path": './MSEED/',
    "SGY_path": './SGY/',
    
    "experiment_path": './DAS-DATA/',

    "events_path": './EVENTS/',
    
    "mseed_land_stations": 'land-stations.mseed',
    "figures_path": './FIGURES/'
}

clock_offset=150

# default parameters for frequency analysis
freq_params={
    "decimation_factor": [8], #[5]
    "land_station_min_freq": 1.,
    "land_station_max_freq": 50.,
}

station_coordinates = {
    "PMAR": {"latitude": 32.72, "longitude": -16.91, "channels": "HH*"},
    "PMPST": {"latitude": 33.08, "longitude": -16.33, "channels": "HH*"},
    "PMPS": {"latitude": 33.06, "longitude": -16.33, "channels": "HN*"},
    "PMOZ": {"latitude": 32.82, "longitude": -17.20, "channels": "HN*"},
    "FUL": {"latitude": 32.65, "longitude": -16.89, "channels": "SH*"},
    "FIMU":{"latitude": 32.29856, "longitude": -16.62694},
}

# default acquisition parameters for OptoDAS
acquisition_params={
    "roi_step": 5,
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
                "net_no_stations": 11293,
                "sta_name_prefix": "GeoLAB fibre. DAS",
                "inventory_source": "maloureiro@ciencias.ulisboa.pt",
                "station_prefix": "G",
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
                "net_no_stations": 11293,
                "sta_name_prefix": "GeoLAB fibre. DAS",
                "inventory_source": "maloureiro@ciencias.ulisboa.pt",
                "station_prefix": "G",
             }
}

channel_ranges={
    "seismic": {"channel_min":  [450] ,
                "channel_max":  [11293] ,
                "channel_step": [1]} ,

    "seismic5": {"channel_min":  [1000] ,
                 "channel_max":  [11283] ,
                 "channel_step": [5]} ,

    "seismic10": {"channel_min":  [1000] ,
                  "channel_max":  [11283] ,
                  "channel_step": [10]} ,

    "cross_correlation": {"channel_min":  [3000] ,
                          "channel_max":  [11000] ,
                          "channel_step": [5]} ,
    
    "internal_tides": {"channel_min":  [500] ,
                       "channel_max":  [11000] ,
                       "channel_step": [100]} ,
    
    "mseed_A": {"channel_min":  [200] ,
                "channel_max":  [4000] ,
                "channel_step": [20]} ,
    
    "coastal": {"channel_min":  [450] ,
                "channel_max":  [2500] ,
                "channel_step": [1]} ,
    
    "land": {"channel_min":  [0] ,
             "channel_max":  [800] ,
             "channel_step": [1]} ,
    
    "spectrogram": {"channel_min":  [200,1200,2200,3200,4200,5200,6200,7200,8200,9200,10200,11200] , # channel ranges need not be contiguous
                    "channel_max":  [305,1205,2205,3205,4205,5205,6205,7205,8205,9205,10205,11205] ,
                    "channel_step": [10 ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1    ,1    ]} ,
}

offset_ranges={
    "2.3-8km + 10-15km":{
           "offset_min": [2300,10000] ,
           "offset_max": [8000,15000] ,  
          "offset_step": [5   ,5    ]},
}
