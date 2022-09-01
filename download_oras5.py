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
# Download files
##########################################################################################
cdo = Cdo()    # Parameters
dwnld_files = []
for data_params in cfg.data_params_all:
    variable = data_params['variable']
    print(f"Download variable {variable}")

    # Vertical resolution depends on the variable type
    if variable in ['meridional_velocity', 'potential_temperature', 
                    'salinity', 'zonal_velocity']:
        levels = 'all_levels'
    else:
        levels = 'single_level'

    # Select years
    if 'time_range' in data_params:
        starty = data_params['time_range'][0]
        endy = data_params['time_range'][1]
    else:
        starty = 1958
        endy = 2022
    print(f"Download year {starty} to {endy}")
    years = np.arange(starty, endy+1, 1)

    filelist = []
    for year in years:
        dirpath = cfg.lpaths['raw_data_dir'] + f"/oras5/{variable}"

        prefix = (dirpath + f"/{variable}_oras5_{levels}_{year}")
        fname = prefix + ".zip"

        # Create directory
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)


        # Change between consolidated (..-2014) and operational (2015-..) product type
        if year < 2015:
            product_type = 'consolidated'
        else:
            product_type = 'operational'


        if (cfg.overwrite) | (not os.path.exists(prefix + ".nc")):
            c = cdsapi.Client()

            c.retrieve(
                'reanalysis-oras5',
                {
                    'product_type': product_type,
                    'format': 'zip',
                    'vertical_resolution': levels,
                    'variable': variable,
                    'year': str(year),
                    'month': [
                        '01', '02', '03',
                        '04', '05', '06',
                        '07', '08', '09',
                        '10', '11', '12',
                    ],
                },
                fname
            )

            # unzip
            print("Unzip files.")
            with ZipFile(fname) as z:
                print(f"Make Dir: {prefix}")
                os.makedirs(prefix)
                z.extractall(prefix)

            month_files=[]
            for f in os.scandir(prefix):
                month_files.append(f.path)

            print("Unzip files.")
            cdo.mergetime(options='-b F32 -f nc', input=month_files,
                          output=prefix + ".nc")

            del c

        filelist.append(prefix + ".nc")

    # Merge yearly files
    prefix_yrange = (dirpath + f"/{variable}_oras5_{levels}_{starty}_{endy}")
    fname_yrange = prefix_yrange + "_raw.nc"
    cdo.mergetime(options='-b F32 -f nc', input=filelist,
                  output=fname_yrange)
    dwnld_files.append(dict(
        fname=fname_yrange,
        prefix=prefix_yrange,
        variable=variable,
    ))

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