import sys, os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as ctp

import geo_preprocessing as gp

def store_monthly_anomalies():
    """Calculate and store monthly anomalies from mon mean."""
    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname +"era5_sst_2000-2019_mon_mean.nc"
    ano_fname = dirname + "era5_sst_2000-2019_mon_anomalies.nc"

    ds = xr.open_dataset(fname)
    # reduce map resolution
    ds_prepro = gp.reduce_map_resolution(ds, lon_factor=20, lat_factor=20)
    # obtain anomalies
    ds_prepro = gp.get_monthly_anomalies(ds_prepro)
    gp.save_to_file(ds_prepro, ano_fname)


    ds = xr.open_dataset(ano_fname)
    # change coordinates
    ds_prepro = gp.set_antimeridian2zero(ds)
    anomalies = ds_prepro['analysed_sst']

    # plotting
    dmap = anomalies[0,:,:]
    gp.plot_map(dmap, central_longitude=180)

    return

def get_nino34_anomaly_lit():
    # Read nino34 from lit
    nino34_file = open("/home/jakob/climate_data/local/climate_indices/nino34_noaa.dat")

    buff = []
    for line in nino34_file.readlines()[2:]:
        buff.append(line.split())
    file_arr = np.array(buff)

    years = file_arr[:, 0]
    months = np.array(file_arr[:, 1], dtype=int)
    nino34_ano = np.array(file_arr[:, 9], dtype=float)

    time = []
    for i, y in enumerate(years):
        time.append(np.datetime64('{}-{:02d}'.format(y, months[i]), 'D'))
    time = np.array(time)

    da = xr.DataArray(data=nino34_ano, name='nino3.4 lit', coords={"time": time}, dims=["time"])

    return da

if __name__ == "__main__":
    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname +"era5_sst_2000-2019_mon_anomalies.nc"

    ds = xr.open_dataset(fname)
    # change coordinates
    ds_prepro = gp.set_antimeridian2zero(ds)
    anomalies = ds_prepro['analysed_sst']

    # average time series
    nino34_mon_mean, nino34_mon_var = gp.get_time_series_area(anomalies, 
        lon_range=gp.get_antimeridian_coord([-170, -120]),
        lat_range=[-5, 5]
    ) 

    # rolling mean over 5 months
    nino34 = nino34_mon_mean.rolling(time=5, center=True).mean()# .dropna("time")

    # plotting monthly mean with std
    plt.figure()
    plt.plot(nino34_mon_mean.time, nino34_mon_mean.data, 
        'o--', color='b',
        label='mon. mean Nino3.4'
    )
    plt.fill_between(nino34_mon_var.time,
        nino34_mon_mean.data - nino34_mon_var.data,
        nino34_mon_mean.data + nino34_mon_var.data,
        alpha=0.3, edgecolor='blue', facecolor='blue',
        label='Variance'
    )

    plt.plot(nino34.time, nino34.data, 
        '-', color='k',
        label='Nino3.4 calc'
    )

    # literature nino3.4
    nino34_lit = get_nino34_anomaly_lit().sel(time=slice("2000-01", "2019-12"))
    plt.plot(nino34_lit.time, nino34_lit.data,
        '-', color='r',
        label='Nino3.4 lit'
    )

    plt.legend(loc=0)
    plt.show()


    
        