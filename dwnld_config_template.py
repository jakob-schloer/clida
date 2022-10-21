#!/usr/bin/env python
import os
PATH = os.path.dirname(os.path.abspath(__file__))

# local paths (imput your own)
lpaths = {
	"raw_data_dir": PATH + "/climate_data"
}
overwrite = True


# CMIP6 specification
# ======================================================================================
data_params_all = [
    {'experiment_id': 'historical', 'source_id': 'CESM2',
     'member_id': 'r1i1p1f1',  # TODO: remove to download all members
     'table_id': 'Amon', 'variable_id': 'ts'},
    {'experiment_id': 'historical', 'source_id': 'CESM2',
     'member_id': 'r1i1p1f1',  # TODO: remove to download all members
     'table_id': 'fx', 'variable_id': 'tas'},
]

pp_params_all = [{
    'lon_range':180,    # 180 or 360 - do you want your output file to count lon -180:180 or 0:360?
    'lon_origin':0,
#    'time': {'historical': ['1850-01-01', '2014-12-31']}, # specify a time range for an experiment_id
#    'lat':[-90, 90],   # Cut lat range
#    'lon':[-180, 180], # Cut lon range
#    'fn_suffix':'',    # added to end of filename when saving
}]  # set origin (first lon value) of pre-processed grid.


# ERA5 
# ======================================================================================
data_params_all = [
    {'variable': 'sea_surface_temperature', 
     'resolution': 'hourly', # 'hourly', 'monthly'
     'plevels': [1000], 
     'time_range': [2017, 2019]}
]
pp_params_all = [{
    'grid_step' : 1,
    'lat': [-60, 60], 
    'lon': [-150, 150], 
    'time_average': 'day', 
}]

# ORAS5 
# ======================================================================================
data_params_all = [
    {'variable': 'sea_surface_temperature', 'time_range': [1958, 2022]}
    {'variable': 'meridional_velocity', 'time_range': [1958, 2022]}
    {'variable': 'zonal_velocity', 'time_range': [1958, 2022]}
]
pp_params_all = [{
    'grid_step' : 1,
    'lat': [-60, 60], 
    'lon': [-150, 150], 
}]

# ERSSTv5, HadISST, COBE2, 
# ======================================================================================
data_params_all = [
    {'variable': 'sst'},
    {'variable': 'ssh'},
    {'variable': 'ucur'}, 
    {'variable': 'vcur'},
]

pp_params_all = [{
    'grid_step' : 1,
    'lat': [-60, 60], 
    'lon': [-150, 150], 
    'time_average': 'month', 
}]
