"""Animation of SST fields.

    #TODO: not working

"""

import sys, os
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
#matplotlib.use("Agg")

import cartopy.crs as ccrs
import cartopy as ctp

import geo_preprocessing as gp

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


dirname = "/home/jakob/climate_data/local/era5/"
ano_red_fname = dirname + "era5_sst_2018-2019_mon_anomalies_redRes.nc"
ds = xr.open_dataset(ano_red_fname)
da = ds['analysed_sst']

fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.coastlines()
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)    
ax.add_feature(ctp.feature.RIVERS)
ax.add_feature(ctp.feature.BORDERS, linestyle=':')

colormap='RdBu'
cmap=plt.get_cmap(colormap)
vmin = -5.0 # np.nanmin(da.data)
vmax = 5.0 # np.nanmax(da.data)

lons = da['lon']
lats = da['lat']

def animate(frame):
    ax.clear()
    ax.coastlines()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)    
    ax.add_feature(ctp.feature.RIVERS)
    ax.add_feature(ctp.feature.BORDERS, linestyle=':')
    dmap = da.isel(time=frame)
    cont = plt.contourf(lons, lats, dmap, 60, cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    return cont

ani = animation.FuncAnimation(fig, animate, frames=np.arange(0, da.shape[0]), blit=False)
ani.save(dirname+'animation.mp4', writer=writer) 

plt.show()
