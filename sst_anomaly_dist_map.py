import sys, os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import geo_preprocessing as gp

if __name__ == "__main__":
    
    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname + "era5_sst_2000-2019_mon_anomalies_redRes.nc"
    ds = xr.open_dataset(fname) 

    # change coordinates
    ds_prepro = gp.set_antimeridian2zero(ds)
    anomalies = ds_prepro['analysed_sst']
    # cut pacific
    lon_range = gp.get_antimeridian_coord([-70, 120])
    pacific_ano = gp.cut_map_area(anomalies, 
        lon_range=lon_range, lat_range=[-30, 30])
    # obtain mean and std 
    buff = pacific_ano.groupby('time.month')
    pacific_mean = buff.mean(dim='time', skipna=True)
    pacific_std = buff.std(dim='time', skipna=True)

    # plotting
    gp.plot_map(pacific_mean.sel(month=1))
    pacific_5quantile = pacific_mean.sel(month=1) - 1.645 * pacific_std.sel(month=1)
    gp.plot_map(pacific_5quantile)
    pacific_95quantile = pacific_mean.sel(month=1) + 1.645 * pacific_std.sel(month=1)
    gp.plot_map(pacific_95quantile)

    plt.show()