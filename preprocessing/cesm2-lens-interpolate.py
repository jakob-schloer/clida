''' Merge CESM2 LENS data into one file and interpolate to a regular grid. 

@Author  :   Jakob Schlör 
@Time    :   2023/08/18 11:29:48
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
import os, argparse
import numpy as np
import xarray as xr
from scipy.interpolate import griddata
from tqdm import tqdm
from joblib import Parallel, delayed

# Run in ipython mode
if False:
    config = dict(
        path="/mnt/qb/goswami/data/cmip6_lens/month/piControl/CESM2/ssh/",
        var='SSH',
        n_cpus = 8,
        grid_step=1.0
    )
else:
    parser = argparse.ArgumentParser()
    parser.add_argument('-path', '--path', default='./',
                        type=str, help='Folder or filepath.')
    parser.add_argument('-var', '--var', default='SSH',
                        type=str, help='Variable name.')
    parser.add_argument('-cpus', '--n_cpus', default=8,
                        type=int, help='Number of cpus.')
    parser.add_argument('-grid', '--grid_step', default=1.0,
                        type=float, help='Grid step.')
    config = vars(parser.parse_args())

var = config['var']
grid_step = config['grid_step']
new_lats = np.arange(-70, 71, grid_step)
new_lons = np.arange(0, 360, grid_step)

# %%
# Load and merge data
# ======================================================================================
print("Load and merge data")
orig_ds = xr.open_mfdataset(os.path.join(config['path'], '*.nc'), chunks={'time': 240},
                            combine='by_coords')
orig_da = orig_ds[var]
orig_lats = orig_da['TLAT'].data
orig_lons = orig_da['TLONG'].data

if var == 'SST':
    orig_da = orig_da.isel(z_t=0)

# %%
# Interpolate to regular grid
# ======================================================================================
def interp_to_regular_grid(da, timeidx, orig_lats, orig_lons,
                           lat_grid, lon_grid, method='linear'):
    """Interpolation of dataarray to a new set of points.

    Args:
        da (xr.Dataarray): Dataarray with coords 
        points_origin (np.ndarray): Array of origin locations.
        points_grid (np.ndarray): Array of locations to interpolate on.
        i (int): Index of time in case of parralelization.

    Returns:
        i (int): Index of time
        values_grid_flat (np.ndarray): Values on new points.
    """
    if timeidx is None:
        orig_data = da.data
    else:
        orig_data = da.isel(time=timeidx).data
    values_grid = griddata(
        (orig_lats.ravel(), orig_lons.ravel()),
        orig_data.ravel(),
        (lat_grid, lon_grid),
        method=method
    )
    return values_grid, timeidx

# Create array of new grid points 
lon_grid, lat_grid = np.meshgrid(new_lons, new_lats)

# Interpolate each time step in parallel
print("Interpolate to regular grid")
n_processes = len(orig_da['time'])
results = Parallel(n_jobs=config['n_cpus'])(
    delayed(interp_to_regular_grid)(
        orig_da, timeidx, orig_lats, orig_lons, lat_grid, lon_grid, method='linear'
    ) for timeidx in tqdm(range(n_processes))
)

# Store interpolated dataarray 
interpolated_data = np.array([r[0] for r in results])
timeids = np.array([r[1] for r in results])
timepoints = orig_da['time'].data[timeids]

interp_da = xr.DataArray(
    data=np.array(interpolated_data),
    dims=['time', 'lat', 'lon'],
    coords=dict(time=timepoints, lat=new_lats, lon=new_lons),
    name=orig_da.name
)

# %%
# Save interpolated data
# ======================================================================================
outpath = os.path.join(config['path'], 
                       f'b.e21.B1850.f09_g17.CMIP6-piControl.001.{var}_interp_gr{grid_step}.nc')
print(f"Save interpolated data to {outpath}")
interp_da.to_netcdf(outpath)
# %%
