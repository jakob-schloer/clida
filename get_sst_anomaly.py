import sys, os
import numpy as np
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from cartopy import config
import cartopy.crs as ccrs
import cartopy as ctp

import geo_preprocessing as gp


def get_monthly_anomalies(ds, output_file=False):
    """Calculate monthly anomalies."""
    climatology = ds.groupby("time.month").mean("time")
    anomalies = ds.groupby("time.month") - climatology

    if output_file is True:
        try:
            anomalies.to_netcdf(output_file)
        except OSError:
            print("Could not write to file!") 

    return anomalies

def get_time_series_area(da, lat_range, lon_range):
    """Get the weighted mean time-series of an area with std."""
    # Cut area of interest
    da_area = da.sel(lat=slice(np.min(lat_range), np.max(lat_range)),
                    lon=slice(np.min(lon_range), np.max(lon_range)))
    
    # Weighted average
    weights = np.cos(np.deg2rad(da_area.lat))
    da_area.weighted()

    
    return


def plot_basic_map(central_longitude=0,  lat_range=[-90,90], lon_range=[-180,180], savepath='./plots/basic_map_plot.pdf' ):
    """Create basic map."""
    fig, ax = plt.subplots(figsize=(9,6))

    # Possible projections are e.g.
    # projections = [ccrs.PlateCarree(), ccrs.Robinson(), ccrs.Mercator(), ccrs.Orthographic(),]

    # Extent
    extent=lon_range+lat_range

    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])
    projection=ccrs.PlateCarree(central_longitude=central_longitude)
    ax=plt.axes(projection=projection)
    ax.coastlines()
    ax.set_global()

    ax.set_extent(extent, ccrs.PlateCarree(central_longitude=central_longitude))

    # Features
    ax.add_feature(ctp.feature.RIVERS)
    ax.add_feature(ctp.feature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    return fig, ax, projection
 

def plot_data_on_global_map(data, lon, lat, lat_range=[-90,90], lon_range=[-180,180], vmin=None, vmax=None, central_longitude=0):
    """Plot colour map on global map"""
    fig, ax, projection = plot_basic_map(central_longitude=central_longitude, lat_range=lat_range, lon_range=lon_range)
    # Data plots
    Lon, Lat = np.meshgrid(lon.values,lat.values)
    x,y=Lon+.25,Lat+.25
    #var_data_mean.plot(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.4})
    colormap='RdBu'
    cmap=plt.get_cmap(colormap)
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)
    step=9
    epsilon=0.001
    ticks=ticks=np.linspace(vmin,vmax,step)
    print("ticks:" , ticks)
    norm = mpl.colors.BoundaryNorm(np.linspace(vmin,vmax+epsilon,2*step), cmap.N)
    cbar=ax.pcolormesh(x,y,data[:], transform=ccrs.PlateCarree(), cmap=cmap,vmin=vmin, vmax=vmax, norm=norm, )

    fig.colorbar(cbar, extend='both', orientation='horizontal', label='SST anomalies [K]', shrink=0.8, ticks=ticks )

    return fig, ax, projection


def plot_animate_maps(da):
    """Create animation of sst map.
        TODO: This does not work yet
        This works only in __main__
    
    """

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    colormap='RdBu'
    cmap=plt.get_cmap(colormap)
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)

    def animate(frame):
        lons = da['lon']
        lats = da['lat']
        data = da.data[frame, :, :]
        cont = plt.contourf(lons, lats, data, 60, cmap=cmap,vmin=vmin, vmax=vmax,
                       transform=ccrs.PlateCarree())
        return cont

    ani = FuncAnimation(fig, animate, frames=np.arange(0, da.shape[0]))
    ani.save('animation.mp4') # TODO: enable animation saving
    return


if __name__ == "__main__":

    dirname = "/home/jakob/climate_data/local/era5/"
    fname = dirname +"era5_sst_2018-2019_mon_mean.nc"
    anom_fname = dirname + "era5_sst_2018-2019_mon_anomalies.nc"
#   anomalies = get_monthly_anomalies(ds, anom_fname)

    ds = xr.open_dataset(anom_fname)
    da_ano = ds['analysed_sst']

    da = da_ano.coarsen(lon=4).mean().coarsen(lat=4).mean()

    data = da.data[0,:,:]
    lats = da.coords['lat']
    lons = da.coords['lon']
    plot_data_on_global_map(data, lons, lats)

    plt.show()


    
        