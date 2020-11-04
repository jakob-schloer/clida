import sys, os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import geo_preprocessing as gp

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

    idx = np.random.randint(0, high=240)
    gp.plot_map(anomalies[-48], central_longitude=180, color='RdBu_r')
    plt.show()