''' Download ORAS5 reanalysis dataset.

@Author  :   Jakob Schlör
@Time    :   2022/08/12 09:08:23
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
import sys
import os
import numpy as np
import xarray as xr
from cdo import Cdo
from zipfile import ZipFile
import cdsapi
import multiprocessing as mpi

import utils as ut

# Get config file
import dwnld_config as cfg

# %%
# Read files
dwnld_files = [dict(
    variable='sea_surface_height',
    fname="/home/jakob/mount_volume/data/oras5/sea_surface_height/sea_surface_height_oras5_single_level_2022.nc",
    prefix="climate_data/ssh_oras5_single_level_2022",
)]

# %%
# Preprocess downloaded files
# ======================================================================================
var_spec = {
    'sea_surface_temperature': dict(vname='sosstsst', new_vname='sst'),
    'sea_surface_height': dict(vname='sossheig', new_vname='ssh'),
    'ocean_heat_content_for_the_upper_300m': dict(vname='sohtc300', new_vname='t300'),
}

for i, f_dwnld in enumerate(dwnld_files):
    print(f"Preprocess {f_dwnld['fname']}", flush=True)
    if len(cfg.pp_params_all) > 1:
        raise ValueError(
            "For ORAS5 we only allow all files to be processed equally!")
    else:
        # Process all downloaded files equally
        pp_params = cfg.pp_params_all[0]

    if f_dwnld['variable'] not in list(var_spec.keys()):
        raise ValueError(f"Add variable shortname to variables dictionary.")

    vname = var_spec[f_dwnld['variable']]['vname']
    prefix = f_dwnld['prefix']

    # Open file
    ds = xr.open_dataset(f_dwnld['fname'])
    da = ds[vname]

    # Interpolate tripolar grid on mercato grid
    grid_step = pp_params['grid_step'] if 'grid_step' in list(pp_params.keys()) else 1
    print(f"Interpolate tripolar grid on mercator grid with res {grid_step}",
          flush=True)
    init_lat = np.arange(
        -90, 90, grid_step
    )
    init_lon = np.arange(
        -180, 180, grid_step 
    )
    grid = {'lat': init_lat, 'lon': init_lon}
    try:
        n_cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
        print(f"SLURM_CPUS_PER_TASK: {os.environ['SLURM_CPUS_PER_TASK']}")
    except:
        print("Could not find environment varialble 'SLURM_CPUS_PER_TASK'.")
        n_cpus = mpi.cpu_count()

    da = ut.interp_points2mercato(da, grid=grid, n_cpus=n_cpus)
    prefix += f"_{grid_step}x{grid_step}"

    # Other preprocessing
    da = ut.check_dimensions(da)

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

    # Save to file
    new_vname = var_spec[f_dwnld['variable']]['new_vname']
    ut.save_to_file(da, prefix + ".nc", var_name=new_vname)
# %%
