import sys, os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import geo_preprocessing as gp

def plot_monthly_maps(maps, title, vmin, vmax, cmap):

    fig = plt.figure(figsize=(20,10))
    for month in np.arange(1,13,1):
        ax = plt.subplot(4,3, month, projection=ccrs.PlateCarree(central_longitude=180))
        ax = gp.plot_map(maps.sel(month=month), central_longitude=180,
                          vmin=vmin, vmax=vmax, ax=ax, color=cmap)
        ax.set_title('month: {}'.format(month))
    
    plt.suptitle(title)


if __name__ == "__main__":
    
    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname + "era5_sst_2000-2019_mon_anomalies.nc"
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
    pacific_5quant = buff.quantile(0.05, dim='time', skipna=True)
    pacific_95quant = buff.quantile(0.95, dim='time', skipna=True)

    # plotting
    vmin = np.nanmin(pacific_std.data)
    vmax = np.nanmax(pacific_std.data)
    plot_monthly_maps(pacific_std, 'Standard deviation of SST-anomalies', vmin, vmax, cmap='YlGn')
    
    vmin = np.nanmin(pacific_5quant.data)
    vmax = np.nanmax(pacific_95quant.data)
    plot_monthly_maps(pacific_5quant, '5% quantile of SST-anomalies', vmin, vmax, cmap='RdBu_r')
    plot_monthly_maps(pacific_95quant, '95% quantile of SST-anomalies', vmin, vmax, cmap='RdBu_r')

    
    plt.show()