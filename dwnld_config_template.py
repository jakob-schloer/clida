#!/usr/bin/env python
import os
PATH = os.path.dirname(os.path.abspath(__file__))

# local paths (imput your own)
lpaths = {
	"raw_data_dir": PATH + "/climate_data"
}

# CMIP6 specification
data_params_all = [
    {'experiment_id': 'historical', 'source_id': 'CESM2',
     'member_id': 'r1i1p1f1',  # TODO: remove to download all members
     'table_id': 'Amon', 'variable_id': 'ts'},
]

# Set parameters for spatiotemporal preprocessing
pp_params_all = [{
    'time': {
        'historical': ['1850-01-01', '2014-12-31'],
        'piControl': ['0001-01-01', '1200-12-31']},
    },                  # make sure to specify the experiment id for each time range
    'lon_range':180,    # 180 or 360 - do you want your output file to count lon -180:180 or 0:360?
    'lon_origin':0,
#    'lat':[-90, 90],   # Cut lat range
#    'lon':[-180, 180], # Cut lon range
#    'fn_suffix':'',    # added to end of filename when saving
}]  # set origin (first lon value) of pre-processed grid.
