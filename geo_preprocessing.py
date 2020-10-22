"""Collection of functions to preprocess climate data."""

import sys, os
import numpy as np
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from cartopy import config
import cartopy.crs as ccrs
import cartopy as ctp

def save_to_file(ds, filepath):
    """Save dataset or dataarray to file."""
    if os.path.exists(filepath):
        print("File" + filepath + " already exists!")
    else:
        try:
            ds.to_netcdf(filepath)
        except OSError:
            print("Could not write to file!") 
    return


def get_monthly_anomalies(ds):
    """Calculate monthly anomalies."""
    climatology = ds.groupby("time.month").mean("time")
    anomalies = ds.groupby("time.month") - climatology
    print('Created monthly anomalies!')
    return anomalies


def get_antimeridian_coord(lons):
    """Change of coordinates from normal to antimeridian."""
    lons = np.array(lons)
    lons_new = np.where(lons < 0, (lons % 180),(lons % 180 - 180)) 
    return lons_new

def set_antimeridian2zero(ds):
    """Set the antimeridian to zero.

    Easier to work with the pacific then.
    """
    # Roll data such that it is centered around the antimeridian
    ds_rolled = ds.roll(lon=(ds.dims['lon'] // 2))
    # Change lon coordinates
    lons = ds_rolled.lon
    lons_new = get_antimeridian_coord(lons)
    ds_rolled = ds_rolled.assign_coords(
        lon=lons_new
    )
    print('Set the antimeridian to the new longitude zero.')
    return ds_rolled 

def reduce_map_resolution(ds, lon_factor, lat_factor):
    """Reduce resolution of map of dataarray."""
    ds_coursen = ds.coarsen(lon=lon_factor).mean().coarsen(lat=lat_factor).mean()
    print('Reduced the resolution of the map!')
    return ds_coursen


def cut_map_area(ds, lon_range, lat_range):
    """Cut an area in the map."""
    ds_area = ds.sel(
        lon=slice(np.min(lon_range), np.max(lon_range)),
        lat=slice(np.min(lat_range), np.max(lat_range))
    )
    return ds_area


def get_time_series_area(da, lon_range, lat_range):
    """Get the weighted mean time-series of an area with std.
    
        The weights are put on the latitudes.
        # TODO: do I have to take care of NaNs?
    """
    # Cut area of interest
    da_area = cut_map_area(da, lon_range, lat_range) 

    # weights of latitudes
    weights = np.cos(np.deg2rad(da_area.lat))

    # weighted mean
    weight_arr = da.weighted(weights)     
    weighted_mean = weight_arr.mean(('lon', 'lat'), skipna=True)
    # weighted variance
    buff = (da - weighted_mean)**2
    weight_arr = buff.weighted(weights)
    weighted_var = weight_arr.mean(('lon', 'lat'), skipna=True)
    
    return weighted_mean, weighted_var


##############################################################################
# Plotting routines from here on
##############################################################################
def plot_map(dmap, central_longitude=0):
    """Simple map plotting using xArray."""
    ax = plt.subplot(projection=ccrs.PlateCarree(central_longitude=central_longitude))

    dmap.plot.pcolormesh("lon", "lat", ax=ax, infer_intervals=True)
    ax.coastlines()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)    
    ax.add_feature(ctp.feature.RIVERS)
    ax.add_feature(ctp.feature.BORDERS, linestyle=':')

    return ax


def plot_animate_maps(da):
    """Create animation of sst map.
        TODO: This does not work yet
        This works only in __main__
    
    """

    fig = plt.figure()
    ax = plt.subplot(projection=ccrs.PlateCarree(central_longitude=central_longitude))
    ax.coastlines()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)    
    ax.add_feature(ctp.feature.RIVERS)
    ax.add_feature(ctp.feature.BORDERS, linestyle=':')

    def animate(frame):
        lons = da['lon']
        lats = da['lat']
        dmap = da.isel(time=frame)
#        cont = plt.contourf(lons, lats, dmap, 60, cmap=cmap,vmin=vmin, vmax=vmax,
#                       transform=ccrs.PlateCarree())
        cont = dmap.plot.pcolormesh("lon", "lat", ax=ax, infer_intervals=True)
        return cont

    ani = FuncAnimation(fig, animate, frames=np.arange(0, da.shape[0]))
    ani.save('animation.mp4') # TODO: enable animation saving
    return