import sys, os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as ctp

import geo_preprocessing as gp


if __name__ == "__main__":

    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname +"era5_sst_2000-2019_mon_mean.nc"
    ano_fname = dirname + "era5_sst_2000-2019_mon_anomalies.nc"
    ano_red_fname = dirname + "era5_sst_2000-2019_mon_anomalies_redRes.nc"

    ds = xr.open_dataset(fname)
    # obtain anomalies
    ds_prepro = gp.get_monthly_anomalies(ds)
    gp.save_to_file(ds_prepro, ano_fname)
    # reduce map resolution
    ds_prepro = gp.reduce_map_resolution(ds_prepro, lon_factor=20, lat_factor=20)
    # save to file
    gp.save_to_file(ds_prepro, ano_red_fname)


    ds = xr.open_dataset(ano_red_fname)
    # change coordinates
    ds_prepro = gp.set_antimeridian2zero(ds)
    anomalies = ds_prepro['analysed_sst']

    # average time series
    nino35_mean, nino35_var = gp.get_time_series_area(anomalies, 
        lon_range=gp.get_antimeridian_coord([-170, -120]),
        lat_range=[-5, 5]
    ) 

    plt.plot(nino35_mean.time, nino35_mean.data, 
        'o-', color='k',
        label='Nino3.5'
    )
    plt.fill_between(nino35_var.time,
        nino35_mean.data - nino35_var.data,
        nino35_mean.data + nino35_var.data,
        alpha=0.3, edgecolor='blue', facecolor='blue',
        label='Variance'
    )
    plt.legend(loc=0)

    # plotting
#    dmap = anomalies[0,:,:]
#    gp.plot_map(dmap, central_longitude=180)


    plt.show()


    
        