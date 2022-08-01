'''Preprocess CMIP6 data

Script for downloading and saving CMIP6 files, with ability to subset by time and space. CMIP6 data is lazily loaded directly from the cloud, using the Pangeo - Google Cloud Public Dataset Program collaboration (more info [here](https://medium.com/pangeo/cmip6-in-the-cloud-five-ways-96b177abe396)).
For each model, files are placed in a subdirectory of the `raw_data_dir` set in the `dwnld_config.py` file `[raw_data_dir]/[model_name]/`. If this subdirectory doesn't yet exist, it is created.

Adopted from Kevin Schwarzwald.

Script is based on Pangeo documentation:
https://pangeo-data.github.io/pangeo-cmip6-cloud/accessing_data.html#opening-a-single-zarr-data-store

@Author  :   Jakob Schlör 
@Time    :   2022/07/27 10:32:30
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
# Packages
##########################################################################################
import xarray as xr
import pandas as pd
import numpy as np
import cftime
import re
from operator import itemgetter  # For list subsetting but this is idiotic
import gcsfs
import os
import warnings


def fix_lons(ds, subset_params):
    """
    This function fixes a few issues that show up when dealing with 
    longitude values. 

    Input: an xarray dataset, with a longitude dimension called "lon"

    Changes: 
    - The dataset is re-indexed to -180:180 or 0:360 longitude format, 
      depending on the subset_params['lon_range'] parameter
    - the origin (the first longitude value) is changed to the closest 
      lon value to subset_params['lon_origin'], if using a 0:360 range. 
      In other words, the range becomes [lon_origin:360 0:lon_origin]. 
      This is to make sure the subsetting occurs in the 'right' direction, 
      with the longitude indices increasing consecutively (this is to ensure
      that subsetting to, say, [45, 275] doesn't subset to [275, 45] or vice-
      versa). Set lon_origin to a longitude value lower than your first subset 
      value.
    """

    if subset_params['lon_range'] == 180:
        # Switch to -180:180 longitude if necessary
        if any(ds.lon > 180):
            ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))
        # Change origin to half the world over, to allow for the
        # longitude indexing to cross the prime meridian, but only
        # if the first lon isn't around -180 (using 5deg as an approx
        # biggest grid spacing). This is intended to move [0:180 -180:0]
        # to [-180:0:180].
        if ds.lon[0] > -175:
            ds = ds.roll(lon=(ds.sizes['lon'] // 2), roll_coords=True)
    elif subset_params['lon_range'] == 360:
        # Switch to 0:360 longitude if necessary
        ds = ds.assign_coords(lon=ds.lon % 360)
        # Change origin to the lon_origin
        ds = ds.roll(
            lon=-((ds.lon // subset_params['lon_origin']) == 1).values.nonzero()[0][0], roll_coords=True)
    return ds


# %%
# Parameters and Variables
##########################################################################################
# Get config file
import dwnld_config as cfg

# Set whether to regrid 360-day calendars to 365-day calendars
regrid_360 = False
overwrite = True

# %%
# Prepare the full query for all the datasets
##########################################################################################

source_calls = np.zeros(len(cfg.data_params_all[0].keys()))

for key in cfg.data_params_all[0].keys():
    if len(np.unique([x[key] for x in cfg.data_params_all])) == 1:
        source_calls[list(cfg.data_params_all[0].keys()).index(key)] = 1

# First get all the ones with the same value for each key
subset_query = ' and '.join([k+" == '"+cfg.data_params_all[0][k]+"'" for k in itemgetter(
    *source_calls.nonzero()[0])(list(cfg.data_params_all[0].keys()))])

# Now add all that are different between subset params - i.e. those that need an OR statement
# These have to be in two statements, because if there's only one OR'ed statement, then the
# for k in statement goes through the letters instead of the keys.
if len((source_calls-1).nonzero()[0]) == 1:
    subset_query = subset_query+' and ('+') and ('.join([' or '.join([k+" == '"+data_params[k]+"'" for data_params in cfg.data_params_all])
                                                         for k in [itemgetter(*(source_calls-1).nonzero()[0])(list(cfg.data_params_all[0].keys()))]])+')'
elif len((source_calls-1).nonzero()[0]) > 1:
    subset_query = subset_query+' and ('+') and ('.join([' or '.join([k+" == '"+data_params[k]+"'" for data_params in cfg.data_params_all])
                                                         for k in itemgetter(*(source_calls-1).nonzero()[0])(list(cfg.data_params_all[0].keys()))])+')'


# %%
# Access google cloud storage links
##########################################################################################
fs = gcsfs.GCSFileSystem(token='anon', access='read_only')
# Get info about CMIP6 datasets
cmip6_datasets = pd.read_csv(
    'https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv'
)
# Get subset based on the data params above (for all search parameters)
cmip6_sub = cmip6_datasets.query(subset_query)

if len(cmip6_sub) == 0:
    warnings.warn('Query unsuccessful, no files found!' +
                  'Check to make sure your table_id matches the domain - for example,' +
                  'SSTs are listed as "Oday" instead of "day"')


# %%
# Process by variable and dataset in the subset
##########################################################################################
for data_params in cfg.data_params_all:
    # Get subset based on the data params above, now just for this one variable
    cmip6_sub = cmip6_datasets.query(
        ' and '.join([k+" == '"+data_params[k] +
                     "'" for k in data_params.keys() if k != 'other'])
    )

    for i, cmip6_sub_row in cmip6_sub.iterrows():
        url = cmip6_sub_row['zstore']

        # Set output filenames
        output_fns = [None]*len(cfg.pp_params_all)
        path_exists = [None]*len(cfg.pp_params_all)

        for subset_params in cfg.pp_params_all:
            # Foldername
            filedir = (cfg.lpaths['raw_data_dir'] + '/'
                       + cmip6_sub_row['table_id'] + '/'
                       + cmip6_sub_row['experiment_id'] + '/'
                       + cmip6_sub_row['source_id'] + '/'
                       + cmip6_sub_row['variable_id'] + '/'
                       )

            # Filename 
            fname = (filedir +
                     cmip6_sub_row['variable_id']+'_' +
                     cmip6_sub_row['table_id']+'_' +
                     cmip6_sub_row['source_id'] + '_' +
                     cmip6_sub_row['experiment_id']+'_' +
                     cmip6_sub_row['member_id']+'_')

            if 'time' in subset_params:
                fname += '-'.join(
                    [re.sub('-', '', t) for t in subset_params['time'][cmip6_sub_row['experiment_id']] ]
                )
            if 'fn_suffix' in subset_params:
                fname += subset_params['fn_suffix']
            fname +='.nc'
            
            output_fns[cfg.pp_params_all.index(subset_params)] = fname

            # Check if path exists
            path_exists[cfg.pp_params_all.index(subset_params)] = os.path.exists(
                output_fns[cfg.pp_params_all.index(subset_params)]
            )

        # Makes overwriting consistent
        if (not overwrite) & all(path_exists):
            warnings.warn('All files already created for ' +
                          cmip6_sub_row['variable_id']+' ' +
                          cmip6_sub_row['table_id']+' ' +
                          cmip6_sub_row['source_id']+' ' +
                          cmip6_sub_row['experiment_id']+' ' +
                          cmip6_sub_row['member_id']+', skipped.')
            continue
        elif any(path_exists):
            if overwrite:
                for subset_params in cfg.pp_params_all:
                    if path_exists[cfg.pp_params_all.index(subset_params)]:
                        os.remove(
                            output_fns[cfg.pp_params_all.index(subset_params)])
                        warnings.warn(
                            'All files already exist for ' +
                            cmip6_sub_row['variable_id']+' ' +
                            cmip6_sub_row['table_id']+' ' +
                            cmip6_sub_row['source_id']+' ' +
                            cmip6_sub_row['experiment_id']+' ' +
                            cmip6_sub_row['member_id'] +
                            ', because OVERWRITE=TRUE theses files have been deleted.'
                        )

        # Open dataset
        print('Download '+ output_fns[cfg.pp_params_all.index(subset_params)])
        ds = xr.open_zarr(fs.get_mapper(url), consolidated=True)

        # Preprocessing of downloaded files
        #################################################################################
        # Rename to lat / lon
        try:
            ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
        except:
            pass

        # same with 'nav_lat' and 'nav_lon' ???
        try:
            ds = ds.rename({'nav_lon': 'lon', 'nav_lat': 'lat'})
        except:
            pass

        # Fix coordinate doubling (this was an issue in NorCPM1,
        # where thankfully the values of the variables were nans,
        # though I still don't know how this happened - some lat
        # values were doubled within floating point errors)
        if 'lat' in ds[data_params['variable_id']].dims:
            if len(np.unique(np.round(ds.lat.values, 10))) != ds.dims['lat']:
                ds = ds.isel(lat=(~np.isnan(ds.isel(lon=1, time=1)[
                             data_params['variable_id']].values)).nonzero()[0], drop=True)
                warnings.warn(
                    'Model ' + ds.source_id +
                    ' has duplicate lat values; attempting to compensate by dropping lat' +
                    ' values that are nan in the main variable in the first timestep'
                )
            if len(np.unique(np.round(ds.lon.values, 10))) != ds.dims['lon']:
                ds = ds.isel(lon=(~np.isnan(ds.isel(lat=1, time=1)[
                             data_params['variable_id']].values)).nonzero()[0], drop=True)
                warnings.warn(
                    'Model '+ds.source_id+' has duplicate lon values; ' +
                    'attempting to compensate by dropping lon values that are nan ' +
                    'in the main variable in the first timestep'
                )

        # Sort by time, if not sorted (this happened with
        # a model; keeping a warning, cuz this seems weird)
        if (ds.time.values != np.sort(ds.time)).any():
            warnings.warn('Model '+ds.source_id +
                          ' has an unsorted time dimension.')
            ds = ds.sortby('time')

        # If 360-day calendar, regrid to 365-day calendar
        if regrid_360 and cmip6_sub_row['table_id'] == 'day':
            if ds.dims['dayofyear'] == 360:
                # Have to put in the compute() because these
                # are by default dask arrays, chunked along
                # the time dimension, and can't interpolate
                # across dask chunks...
                ds = ds.compute().interp(dayofyear=(np.arange(1, 366)/365)*360)
                # And reset it to 1:365 indexing on day of year
                ds['dayofyear'] = np.arange(1, 366)
                # Throw in a warning, too, why not
                warnings.warn('Model ' + ds.source_id +
                              ' has a 360-day calendar; daily values were ' +
                              'interpolated to a 365-day calendar')

        # Save by the subsets desired in cfg.pp_params_all above
        #################################################################################
        for subset_params in cfg.pp_params_all:
            # Make sure this file hasn't already been processed
            if (not overwrite) & path_exists[cfg.pp_params_all.index(subset_params)]:
                warnings.warn(output_fns[cfg.pp_params_all.index(
                    subset_params)]+' already exists; skipped.')
                continue

            # Make sure the target directory exists
            if not os.path.exists(filedir):
                os.makedirs(filedir)
                warnings.warn('Directory ' + filedir+' created!')

            # Fix longitude (by setting it to either [-180:180]
            # or [0:360] as determined by subset_params, and
            # to roll them so the correct range is consecutive
            # in lon (so if you're looking at the Equatorial
            # Pacific, make it 0:360, with the first lon value
            # at 45E).
            try:
                ds_tmp = fix_lons(ds, subset_params)
            except:
                ds_tmp = ds
                warnings.warn(
                    'fix_lons did not work because of the multi-dimensional index'
                )

            # Subset by time as set in subset_params
            time_range = subset_params['time'][cmip6_sub_row['experiment_id']]
            if (ds.time.max().dt.day == 30) | (type(ds.time.values[0]) == cftime._cftime.Datetime360Day):
                # (If it's a 360-day calendar, then subsetting to "12-31"
                # will throw an error; this switches that call to "12-30")
                # Also checking explicitly for 360day calendar; some monthly
                # data is still shown as 360-day even when it's monthly, and will
                # fail on date ranges with date 31 in a month
                ds_tmp = (ds_tmp.sel(time=slice(
                    time_range[0], re.sub('-31', '-30', time_range[1])
                    )))
            else:
                ds_tmp = (ds_tmp.sel(time=slice(*time_range)))

            # Save as NetCDF file
            ds_tmp.to_netcdf(
                output_fns[cfg.pp_params_all.index(subset_params)])

            # Status update
            print(output_fns[cfg.pp_params_all.index(
                subset_params)]+' processed!')

        del ds, ds_tmp, subset_params


# %%
