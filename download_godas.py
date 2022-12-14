''' Download GODAS reanalysis data.

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
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='pottmp.',
        start='1980',
        end='present',
        time_res='month',
        vname='pottmp',
        zlevel=0,
    ),
    'ssh': dict(
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='sshg.',
        start='1980',
        end='present',
        time_res='month',
        vname='sshg',
    ),
    'ucur': dict(
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='ucur.',
        start='1980',
        end='present',
        time_res='month',
        vname='ucur',
        zlevel=0,
    ),
    'vcur': dict(
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='vcur.',
        start='1980',
        end='present',
        time_res='month',
        vname='vcur',
        zlevel=0,
    ),
    'taux': dict(
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='uflx.',
        start='1980',
        end='present',
        time_res='month',
        vname='uflx',
    ),
    'tauy': dict(
        url='https://downloads.psl.noaa.gov/Datasets/godas/',
        prefix='vflx.',
        start='1980',
        end='present',
        time_res='month',
        vname='vflx',
    ),
}
# %%
# Download files
##########################################################################################
dwnld_files = []
for data_params in cfg.data_params_all:
    if data_params['variable'] not in list(variables.keys()):
        warnings.warn(f"Variable '{data_params['variable']}' does not exist!")
        continue

    varspec = variables[data_params['variable']]
 
    # Filepath
    dirpath = (cfg.lpaths['raw_data_dir']
               + f"/GODAS/{varspec['time_res']}/{data_params['variable']}")

    # Create directory
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    # Download data
    urlpath = request.urlopen(varspec['url'])
    string = urlpath.read().decode('utf-8')
    pattern = re.compile(varspec['prefix']+ '[0-9]{4}.nc') 
    remotefilelist = pattern.findall(string)[::2]


    localfilelist = []
    for i, remotefname in enumerate(remotefilelist):
        fname = dirpath + "/" + remotefname
        if (cfg.overwrite) | (not os.path.exists(fname)):
            print(f'Download {remotefname}:', flush=True)
            response = request.urlretrieve(varspec['url'] + remotefname,
                                           fname)
        else:
            print(f"File {fname} exists and will not be overwritten!", flush=True)
        
        localfilelist.append(fname)
        
    # Merge and preprocess files
    filelist = []
    for i, fname in enumerate(localfilelist):
        # Open file
        ds = xr.open_dataset(fname)
        da = ds[varspec['vname']]
        if 'zlevel' in varspec.keys():
            da = da.isel(level=varspec['zlevel'])
        da = ut.check_dimensions(da)
        filelist.append(da)
    da_merge = xr.concat(filelist, dim='time')
    prefix = (dirpath + f"/{data_params['variable']}_godas_{varspec['time_res']}"
              + f"_{varspec['start']}-{varspec['end']}")
    ut.save_to_file(da_merge, prefix + "_raw.nc", var_name=data_params['variable'])

    dwnld_files.append(dict(
        fname=prefix + "_raw.nc",
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
