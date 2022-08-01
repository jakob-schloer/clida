import os,sys
import json
import warnings

import xarray as xr
import xesmf as xe
import h5py as h5
import numpy as np
from tqdm.auto import tqdm
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

dataset_root="../datasets"


# create a model-file dictionary, where values for a given key are the respective model files that need to be rescaled
model_file_dict=dict()
for entry in os.listdir(dataset_root):
     if os.path.isfile(os.path.join(dataset_root, entry)):
         if ".nc" in entry and "pr_" in entry and not "_rescaled" in entry:
             model = entry.split("_")[2]
             if not model in model_file_dict.keys():
                 model_file_dict[model]=[entry]
             else:
                 model_file_dict[model].append(entry)


# The following loop goes through all model files and finds the largest common time-range
t0_max, tn_min = np.datetime64("1970-01-01"), np.datetime64("2020-12-31") #init time range ends
t_range = np.arange(t0_max, tn_min, dtype="datetime64[D]")
print("Finding smallest common time range within the different datasets:")
for model in tqdm(model_file_dict, position=0, desc='Total Progress', leave=False):
    for file in model_file_dict[model]:
        da_ = xr.load_dataset(os.path.join(dataset_root, file))
        # in case dataset time format is cftime.DatetimeNoLeap, convert to numpy datetime format
        if isinstance(da_.indexes['time'], xr.coding.cftimeindex.CFTimeIndex):
            datetimeindex = da_.indexes['time'].to_datetimeindex()
            da_['time'] = datetimeindex
        # find minimum and maximum date-values for given dataset
        t0_ = pd.to_datetime(min(da_.time.values))
        tn_ = pd.to_datetime(max(da_.time.values))
        # update new common time range
        if t0_>=t0_max and tn_<=tn_min:
            # time values are of the form %Y-%m-%dT%H:%M:%s
            # some models dates of the daily values have the hour set to 12, while the rest have the hour set to 0
            # correct the time if hour is set to 12
            if int(t0_.strftime('%H'))==12:
                t_range = da_.time.values - np.timedelta64(12, 'h')
            else:
                t_range = da_.time.values
            t0_max = t0_
            tn_min = tn_

# reformat time range ends to be of the form %Y-%m-%d
t0_max = pd.to_datetime(str(t0_max))
t0_max = t0_max.strftime('%Y-%m-%d')
tn_min = pd.to_datetime(str(tn_min))
tn_min = tn_min.strftime('%Y-%m-%d')

print("The smallest time range found was between {} and {}. All downscaled datasets will span between these dates.\n".format(t0_max, tn_min) )

# Select data inside the calculated time range from the reanalysis dataset
target_file="ECMWF-ERA5_tp_JJA_1979-2021.nc"
da_target=xr.load_dataarray(os.path.join(dataset_root, target_file))
# Save Reanalysis target as h5 file
with h5.File(os.path.join(dataset_root, "rescaled_datasets", "pr_reanalysis_ECMWF-ERA5.h5"), 'w') as f:
     target_data = np.asarray(da_target.sel(time=t_range)*1000) # select time range and turn m/24h into mm/24h
     dset = f.create_dataset("pr_reanalysis", target_data.shape, dtype=np.float64, data=target_data)


print("Downscaling the CMIP6 datasets to the ECMWF-ERA5 Historical Reanalysis dataset...\n")
# init dictionary that lists, given a model name, how many unique instances are available 
available_models=dict()
for model in tqdm(model_file_dict, position=0, desc='Total Progress', leave=True):
    # init model meta data dictionary
    model_meta=dict()
    # get number of unique model instances
    available_models[model]=len(model_file_dict[model])
    for file in tqdm(model_file_dict[model], position=1, desc='Model Progress ({})'.format(model), leave=True):
        filename_base="pr_historical_"+model+"_rescaled"
        model_nr = int(file.split(".")[0].split("_")[-1])

        # check if model HDF5 file already exists and skip any dataset is already part of the database
        if os.path.isfile(os.path.join(dataset_root, "rescaled_datasets", filename_base+".h5")):
            with h5.File(os.path.join(dataset_root, "rescaled_datasets", filename_base+".h5"), 'a') as f:
                if "pr_model_"+str(model_nr) in f.keys(): continue    
        
        # load model dataset
        da_ = xr.load_dataset(os.path.join(dataset_root, file))

        # in case dataset time format is cftime.DatetimeNoLeap, convert to numpy datetime format
        if isinstance(da_.indexes['time'], xr.coding.cftimeindex.CFTimeIndex):
            datetimeindex = da_.indexes['time'].to_datetimeindex()
            da_['time'] = datetimeindex
        # correct hour in date to 0
        t0_ = pd.to_datetime(min(da_.time.values))
        if int(t0_.strftime('%H'))==12:
            da_['time'] = da_['time'] - np.timedelta64(12, 'h')
        # set model attributes
        model_attrs=da_.pr_.attrs
        
        # Regrid model data to target grid
        regridder = xe.Regridder(da_, da_target, "bilinear", periodic=False)
        out_da = regridder(da_)
        # regridded variable name is by defaul "__xarray_dataarray_variable__", set it back to "pr_"
        out_da = out_da.rename({"__xarray_dataarray_variable__":"pr_"})

        # select common time range calculated earlier
        try:
            out_da = out_da.sel(time=t_range)
        except:
            print(da_.time.values)
            print(t_range)
            break
        # convert units from kg m-2 s-1 into mm/24h
        if model_attrs["units"] == "kg m-2 s-1":
            model_data = np.asarray(out_da.pr_)*86400 # turn kg m-2 s-1 into mm/24h
            model_attrs["units"] = "mm/24h"
        else:
            model_data = np.asarray(da_)

        lat_res = abs(np.sum(da_.lat.values[1:] - da_.lat.values[0:-1]))/len(da_.lat.values[0:-1])
        lon_res = abs(np.sum(da_.lon.values[1:] - da_.lon.values[0:-1]))/len(da_.lon.values[0:-1])
        model_attrs["original resolution (lat,lon) [deg]"] = (round(lat_res,2), round(lon_res,2))

        out_da['pr_'] = (('time', 'lat', 'lon'), model_data, model_attrs)

        # save resulting rescaled dataset as a HDF5 dataset and save model attributes to model meta dictionary
        try:
            with h5.File(os.path.join(dataset_root, "rescaled_datasets", filename_base+".h5"), 'a') as f:
                dset = f.create_dataset("pr_model_"+str(model_nr), model_data.shape, dtype=np.float64, data=model_data)
        except ValueError as e:
            print("\n",e)
        model_meta["pr_model_"+str(model_nr)]=model_attrs
    print('\n')
    # save model meta dictionary as a json file
    with open(os.path.join(dataset_root, "rescaled_datasets", filename_base+".json"), 'w') as f:
        f.write(json.dumps(model_meta, default=str, indent=2, separators=(',', ': ')))


# save time, lat and lon coordinates as HDF5 datasets
time_ = [str(t) for t in da_target.time.values]
lat_ = da_target.latitude.values
lon_ = da_target.longitude.values
with h5.File(os.path.join(dataset_root, "rescaled_datasets", "rescaled_coordinates.h5"), 'w') as f:
            dset = f.create_dataset("time", len(time_), data=time_)
            dset = f.create_dataset("lat", len(lat_), dtype=np.float64, data=lat_)
            dset = f.create_dataset("lon", len(lon_), dtype=np.float64, data=lon_)

# save available model dictionary as a json file
with open(os.path.join(dataset_root, "rescaled_datasets", "available_models.json"), 'w') as f:
        f.write(json.dumps(available_models, default=str, indent=2, separators=(',', ': ')))


    
        
        


