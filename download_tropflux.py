''' Unfortunately tropflux data have to be downloaded manually from the following website.

https://incois.gov.in/tropflux/

This script only merges and preprocesses the data.

@Author  :   Jakob Schlör 
@Time    :   2022/07/27 10:32:30
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
# Packages
##########################################################################################
import os
import re
import warnings
import urllib.request as request   
from tqdm import tqdm
import xarray as xr
import pandas as pd
import numpy as np
import utils as ut

# Get config file
import dwnld_config as cfg


# Specification of variables in dataset
variables = {
    'sst': dict(
        folder=cfg.lpaths['raw_data_dir'] + f"/tropflux/month/sst",
        fprefix='sst_tropflux_1m_',
        start=1979,
        end=2018,
        time_res='month',
        vname='sst',
    )
}
# %%
# Download files
##########################################################################################
# Files need to be downloaded manually

# %%
# Merge files
##########################################################################################
dwnld_files = []
for data_params in cfg.data_params_all:

    if data_params['variable'] not in list(variables.keys()):
        warnings.warn(f"Variable '{data_params['variable']}' does not exist!")
        continue
    varspec = variables[data_params['variable']]

    filelist = []
    for i, year in enumerate(np.arange(varspec['start'], varspec['end']+1, 1)):
        fname = varspec['folder'] + f"/{varspec['fprefix']}{year}.nc"
        # Open file
        if os.path.exists(fname):
            ds = xr.open_dataset(fname)
            da = ds[varspec['vname']]
            da = ut.check_dimensions(da)
            filelist.append(da)
        else:
            print(f"File {fname} does not exist.")

    da_merge = xr.concat(filelist, dim='time')
    prefix = (varspec['folder'] + f"/{data_params['variable']}_tropflux_{varspec['time_res']}"
              + f"_{varspec['start']}-{varspec['end']}")
    ut.save_to_file(da_merge, prefix + ".nc", var_name=data_params['variable'])

    dwnld_files.append(dict(
        fname=prefix + ".nc",
        prefix=prefix,
        variable=data_params['variable']
    ))

# %%
# Preprocess downloaded files 
# ======================================================================================
for i, f_dwnld in enumerate(dwnld_files):
    if len(cfg.pp_params_all) > 1:
        # Process files differently
        assert len(cfg.pp_params_all) == len(cfg.data_params_all)
        pp_params = cfg.pp_params_all[i]
    else:
        # Process all downloaded files equally
        pp_params = cfg.pp_params_all[0]
        
    vname = f_dwnld['variable']
    prefix = f_dwnld['prefix']
    # Open file
    ds = xr.open_dataset(f_dwnld['fname'])
    da = ds[vname]
    da = ut.check_dimensions(da)

    # Time averages
    if 'time_average' in list(pp_params.keys()):
        print(f"Resample time by {pp_params['time_average']} and compute mean.", flush=True)
        da = da.resample(time=pp_params['time_average'], label='left').mean()
        da = da.assign_coords(
            dict(time=da['time'].data + np.timedelta64(1, 'D'))
        )
        prefix += f"_tmean_{pp_params['time_average']}"

    # Cut area of interest
    if 'lon' in list(pp_params.keys()):
        lon_range = pp_params['lon'] 
        prefix += f"_lon_{lon_range[0]}-{lon_range[1]}"
    else: 
        lon_range = None
    if 'lat' in list(pp_params.keys()):
        lat_range = pp_params['lat'] 
        prefix += f"_lon_{lat_range[0]}-{lat_range[1]}"
    else: 
        lat_range = None

    if ('lon' in list(pp_params.keys())) or ('lat' in list(pp_params.keys())):
        print(f'Get selected area: lon={lon_range}, lat={lat_range}!', flush=True)
        da = ut.cut_map(
            da, lon_range=lon_range, lat_range=lat_range, shortest=False
        ) 

    # coarse grid if needed
    if 'grid_step' in list(pp_params.keys()):
        print(f"Interpolate grid on res {pp_params['grid_step']}", flush=True)
        da, grid = ut.set_grid(
            da, step_lat=pp_params['grid_step'], step_lon=pp_params['grid_step'],
            lat_range=lat_range, lon_range=lon_range)
        prefix += f"_{pp_params['grid_step']}x{pp_params['grid_step']}"

    # Save to file
    if prefix != f_dwnld['prefix']:
        ut.save_to_file(da, prefix + ".nc", var_name=f_dwnld['variable'])
# %%
