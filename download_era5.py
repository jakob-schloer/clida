''' Download ERA5 reanalysis dataset.

Variable list:
'10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
'2m_temperature','convective_precipitation','convective_rain_rate',
'convective_snowfall','downward_uv_radiation_at_the_surface','evaporation',
'high_cloud_cover','instantaneous_10m_wind_gust','large_scale_precipitation',
'large_scale_rain_rate','large_scale_snowfall','large_scale_snowfall_rate_water_equivalent',
'low_cloud_cover','maximum_2m_temperature_since_previous_post_processing','mean_sea_level_pressure',
'medium_cloud_cover','minimum_2m_temperature_since_previous_post_processing','potential_evaporation',
'precipitation_type','sea_surface_temperature','skin_temperature',
'snow_depth','snow_evaporation','snowfall',
'snowmelt','surface_latent_heat_flux','surface_net_solar_radiation',
'surface_net_thermal_radiation','surface_pressure','surface_sensible_heat_flux',
'surface_solar_radiation_downwards','surface_thermal_radiation_downwards','toa_incident_solar_radiation',
'top_net_solar_radiation','top_net_thermal_radiation','total_cloud_cover',
'total_precipitation','total_sky_direct_solar_radiation_at_surface'


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
    print(f"Download variable {variable}", flush=True)

    # Time resolution
    if 'resolution' in data_params:
        resolution = data_params['resolution']
    else:
        resolution = 'hourly'
    print(f"Download {resolution} data.", flush=True)

    # Pressure levels
    if 'plevels' in data_params:
        plevels = data_params['plevels']
    else:
        # Single pressure level
        plevels = ['sp']

    # Select years
    if 'time_range' in data_params:
        starty = data_params['time_range'][0]
        endy = data_params['time_range'][1]
    else:
        starty = 1940
        endy = 2022
    print(f"Download year {starty} to {endy}", flush=True)
    years = np.arange(starty, endy+1, 1)

    for plevel in plevels:
        print(f"Download pressure level {plevel}", flush=True)
        filelist = []
        for year in years:
            dirpath = cfg.lpaths['raw_data_dir'] + f"/era5" + f"/{resolution}"
            if plevel == 'sp':
                dirpath += f'/single_pressure_level/{variable}/'
            else:
                dirpath += f'/multi_pressure_level/{variable}/{plevel}/'

            prefix = (dirpath + f"/{variable}_era5_{resolution}_{plevel}_{year}")
            fname = prefix + ".nc"

            if (cfg.overwrite) | (not os.path.exists(fname)):
                # Create directory
                if not os.path.exists(dirpath):
                    os.makedirs(dirpath)

                c = cdsapi.Client()

                if resolution == 'hourly':
                    if plevel == 'sp':
                        c.retrieve(
                            'reanalysis-era5-single-levels',
                            {
                                'product_type': 'reanalysis',
                                'format': 'netcdf',
                                'variable': [variable],
                                'year': [str(year)],
                                'month': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                ],
                                'day': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                    '13', '14', '15',
                                    '16', '17', '18',
                                    '19', '20', '21',
                                    '22', '23', '24',
                                    '25', '26', '27',
                                    '28', '29', '30',
                                    '31',
                                ],
                                'time': [
                                    '00:00', '01:00', '02:00',
                                    '03:00', '04:00', '05:00',
                                    '06:00', '07:00', '08:00',
                                    '09:00', '10:00', '11:00',
                                    '12:00', '13:00', '14:00',
                                    '15:00', '16:00', '17:00',
                                    '18:00', '19:00', '20:00',
                                    '21:00', '22:00', '23:00',
                                ],
                            },
                            fname)
                    else:
                        c.retrieve(
                            'reanalysis-era5-pressure-levels',
                            {
                                'product_type': 'reanalysis',
                                'format': 'netcdf',
                                'variable': [variable],
                                'year': [str(year)],
                                'pressure_level': plevel,
                                'month': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                ],
                                'day': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                    '13', '14', '15',
                                    '16', '17', '18',
                                    '19', '20', '21',
                                    '22', '23', '24',
                                    '25', '26', '27',
                                    '28', '29', '30',
                                    '31',
                                ],
                                'time': [
                                    '00:00', '01:00', '02:00',
                                    '03:00', '04:00', '05:00',
                                    '06:00', '07:00', '08:00',
                                    '09:00', '10:00', '11:00',
                                    '12:00', '13:00', '14:00',
                                    '15:00', '16:00', '17:00',
                                    '18:00', '19:00', '20:00',
                                    '21:00', '22:00', '23:00',
                                ],
                            },
                            fname)

                elif resolution == 'monthly':
                    if plevel == 'sp':
                        c.retrieve(
                            'reanalysis-era5-single-levels-monthly-means',
                            {
                                'product_type': 'monthly_averaged_reanalysis',
                                'format': 'netcdf',
                                'variable': [variable],
                                'year': [str(year)],
                                'month': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                ],
                                'time': '00:00',
                            },
                            fname)
                    else:
                        c.retrieve(
                            'reanalysis-era5-pressure-levels-monthly-means',
                            {
                                'product_type': 'monthly_averaged_reanalysis',
                                'format': 'netcdf',
                                'variable': [variable],
                                'year': [str(year)],
                                'pressure_level': plevel,
                                'month': [
                                    '01', '02', '03',
                                    '04', '05', '06',
                                    '07', '08', '09',
                                    '10', '11', '12',
                                ],
                                'time': '00:00',
                            },
                            fname)
                else:
                    ValueError(f"{resolution} does not exist!")
                del c
            else:
                print(f"File {fname} exists and will not be overwritten!", flush=True)

            filelist.append(fname)


        # Merge and preprocess files
        print("Merge files.", flush=True)
        merge_files = []
        for i, fname in enumerate(filelist):
            # Open file
            ds = xr.open_dataset(fname)
            ds = ut.check_dimensions(ds)
            merge_files.append(ds)
        ds_merge = xr.concat(merge_files, dim='time')
        prefix = (dirpath + f"/{variable}_era5_{resolution}_{plevel}"
                  + f"_{starty}-{endy}")
        print(f"Store merged file {prefix}_raw.nc.", flush=True)
        ds_merge.to_netcdf(prefix + "_raw.nc")
            
        dwnld_files.append(dict(
            fname=prefix + "_raw.nc",
            prefix=prefix,
            variable=variable,
        ))

# %%
# Preprocess downloaded files
# ======================================================================================
variables = {
    'sea_surface_temperature': dict(vname='sst'),
    '10m_u_component_of_wind': dict(vname='u10'),
    '10m_v_component_of_wind': dict(vname='v10'),
    'total_precipitation': dict(vname='tp'),
    'geopotential': dict(vname='z'),
    'surface_pressure': dict(vname='sp'),
    '2m_temperature': dict(vname='t2m'),
    'top_net_thermal_radiation' : dict(vname='ttr'),
    'u_component_of_wind': dict(vname='u'),
    'v_component_of_wind': dict(vname='v'),
    'vertical_velocity': dict(vname='w'),
}

for i, f_dwnld in enumerate(dwnld_files):
    print(f"Preprocess {f_dwnld['fname']}", flush=True)
    if len(cfg.pp_params_all) > 1:
        raise ValueError(
            "For ERA5 we only allow all files to be processed equally!")
    else:
        # Process all downloaded files equally
        pp_params = cfg.pp_params_all[0]

    if f_dwnld['variable'] not in list(variables.keys()):
        raise ValueError(f"Add variable shortname to variables dictionary.")
    vname = variables[f_dwnld['variable']]['vname']
    prefix = f_dwnld['prefix']

    # Open file
    ds = xr.open_dataset(f_dwnld['fname'])
    da = ds[vname]
    da = ut.check_dimensions(da)


    # Time averages
    if 'time_average' in list(pp_params.keys()):
        print(f"Resample time by {pp_params['time_average']} and compute mean.",
              flush=True)
        da = da.resample(time=pp_params['time_average'], label='left').mean()
        if pp_params['time_average'] == 'month':
            da = da.assign_coords(
                dict(time=da['time'].data + np.timedelta64(1, 'D'))
            )
        prefix += f"_{pp_params['time_average']}mean"

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
    ut.save_to_file(da, prefix + ".nc", var_name=vname)
