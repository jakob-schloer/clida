import os
import pandas as pd
import numpy as np
import gcsfs
import xarray as xr
import zarr
from tqdm import tqdm
import dask
from dask.diagnostics import ProgressBar

# Script is based on Pangeo documentation: https://pangeo-data.github.io/pangeo-cmip6-cloud/accessing_data.html#opening-a-single-zarr-data-store

dataset_root = "datasets"
df = pd.read_csv('https://cmip6.storage.googleapis.com/pangeo-cmip6.csv')
df.head()

gcs = gcsfs.GCSFileSystem(token='anon', access='read_only')

var_ = ['pr']
for var in var_:
    # for a given var, collect institution IDs (i.e. model names) and the number of available unique models for each institution 
    df_var = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == '"+var+"' & experiment_id == 'historical'")
    df_var['institution_id'].value_counts().to_csv(os.path.join(dataset_root, var+'_historical_data_available_unique_models.csv'))

for var in var_:
    # load df listing model names and the frequency of them for the respective var
    unique_models = pd.read_csv(os.path.join(dataset_root, var+'_historical_data_available_unique_models.csv'))
    table_id = 'day'
    for i_model,model in enumerate(unique_models.iloc[:,0]):
        query_ = "activity_id=='CMIP' & table_id == 'day' & variable_id == '"+var+"' & experiment_id == 'historical' & institution_id =='"+model+"'"
        # for a given model and variable, the following df extracts the relevant entries from df
        df_ta = df.query(query_)
        if df_ta.size==0:
            raise ValueError("query results are empty!")
        for ens in tqdm(range(unique_models.iloc[i_model,1])):
            # for each model instance of a given model name (i.e. NOAA-GFDL) download the respective netCDF file
            # as  {var}_historical_{model name}_ens_{model instance index}.nc
            filename_ = var+'_historical_'+str(model)+'_ens_'+str(ens+1)+'.nc'
            if not os.path.isfile(os.path.join(dataset_root, filename_)):
                # connect to google cloud storage, and then download the respective model files
                zstore = df_ta["zstore"].values[ens] # the google cloud storage url for the respective model instance
        
                # create a mutable-mapping-style interface to the store
                mapper = gcs.get_mapper(zstore)

                # open it using xarray and zarr
                ds_ta = xr.open_zarr(mapper, consolidated=True)

                # calculate latitude and longitude resolution of the given model
                try:
                    lat_values=ds_ta[var].lat.values
                    lon_values=ds_ta[var].lon.values
                    
                except Exception as e:
                    lat_values=ds_ta[var].latitude.values
                    lon_values=ds_ta[var].longitude.values

                lat_res = abs(np.sum(lat_values[1:]-lat_values[:-1]))/len(lat_values[:-1])
                lon_res = abs(np.sum(lon_values[1:]-lon_values[:-1]))/len(lon_values[:-1])
                try:
                    # Select area covering India (compensating for different spatial resolutions), summer season (JJA), years from 1979 to 2021
                    tas_hist_ = ds_ta[var].sel(lat=slice(7.1-lat_res,38.1+lat_res)).sel(lon=slice(67.1-lon_res,98.4+lon_res)).sel(time=ds_ta['time.season']=='JJA').sel(time=slice('1979','2021'))
                    # Set dataset attributes 
                    tas_hist_.attrs = ds_ta[var].attrs
                except Exception as e:
                    print("Error", e)
                    continue
                out_data = tas_hist_
                # Save model data as netCDF file
                with ProgressBar():
                    out_data[var+'_'] = (('time', 'lat', 'lon'), np.asarray(tas_hist_),tas_hist_.attrs)
                    out_data[var+'_'].to_netcdf(os.path.join(dataset_root, filename_))
    
                
                    
                
                