''' Download ORAS5 reanalysis dataset.

@Author  :   Felix Strnad, Jakob Schlör
@Time    :   2022/08/12 09:08:23
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
import sys
import os
from tkinter import Variable
import warnings
import numpy as np
import xarray as xr
from cdo import Cdo
from zipfile import ZipFile
import cdsapi

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


        if cfg.overwrite | os.path.exists(fname):
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

        dwnld_files.append(prefix + ".nc")

    fname_yrange = (dirpath + f"/{variable}_oras5_{levels}_{starty}_{endy}.nc")
    cdo.mergetime(options='-b F32 -f nc', input=dwnld_files,
                  output=fname_yrange)

# %%
# TODO: Preprocessing of oras5. Complication due to tripolar grid
