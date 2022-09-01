"""Collection of functions to preprocess climate data."""
import os
import numpy as np
import scipy as sp
import scipy.interpolate as interpolate
import xarray as xr
from tqdm import tqdm
from joblib import Parallel, delayed


def check_dimensions(ds, sort=True):
    """
    Checks whether the dimensions are the correct ones for xarray!
    """
    dims = list(ds.dims)

    rename_dic = {
        'longitude': 'lon',
        'nav_lon': 'lon',
        'xt_ocean': 'lon',
        'latitude': 'lat',
        'nav_lat': 'lat',
        'yt_ocean': 'lat',
        'time_counter': 'time',
    }
    for c_old, c_new in rename_dic.items():
        if c_old in dims:
            print(f'Rename:{c_old} : {c_new} ')
            ds = ds.rename({c_old: c_new})
            dims = list(ds.dims)

    # Check for dimensions
    clim_dims = ['time', 'lat', 'lon']
    for dim in clim_dims:
        if dim not in dims:
            raise ValueError(
                f"The dimension {dim} not consistent with required dims {clim_dims}!")

    # If lon from 0 to 360 shift to -180 to 180
    if max(ds.lon) > 180:
        print("Shift longitude!")
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))

    if sort:
        print('Sort longitudes and latitudes in ascending order, respectively')
        ds = ds.sortby('lon')
        ds = ds.sortby('lat')

    # if 'time' in ds.dims:
        # ds = ds.transpose('time', 'lat', 'lon')

    return ds


def save_to_file(da, filepath, var_name=None):
    """Save dataset or dataarray to file."""
    if os.path.exists(filepath):
        print("File" + filepath + " already exists!")

    # convert xr.dataArray to xr.dataSet if needed
    if var_name is not None:
        ds = da.to_dataset(name=var_name)
    else:
        ds = da 

    # Store to .nc file
    if os.path.exists(filepath):
        print(f"File '{filepath}' already exists!")
        filepath = filepath + "_new"

    try:
        ds.to_netcdf(filepath)
        print(f"File is stored to '{filepath}'!")
    except OSError:
        print(f"Could not write to '{filepath}'!")

    return None


def time_average(ds, group='1D'):
    """Downsampling of time dimension by averaging.

    Args:
    -----
    ds: xr.dataFrame 
        dataset
    group: str
        time group e.g. '1D' for daily average from hourly data
    """
    ds_average = ds.resample(time=group, label='left').mean(skipna=True)

    # Shift time to first of month
    if group == '1M':
        new_time = ds_average.time.data + np.timedelta64(1, 'D')
        new_coords = {}
        for dim in ds_average.dims:
            new_coords[dim] = ds_average[dim].data
        new_coords['time'] = new_time
        ds_average = ds_average.assign_coords(coords=new_coords)

    return ds_average


def set_grid(ds, step_lat=1, step_lon=1,
             lat_range=None, lon_range=None):
    """Interpolate grid.

    Args:
        ds (xr.Dataset): Dataset or dataarray to interpolate.
            Dataset is only supported for grid_type='mercato'.
        step_lat (float, optional): Latitude grid step. Defaults to 1.
        step_lon (float, optional): Longitude grid step. Defaults to 1.

    Returns:
        da (xr.Dataset): Interpolated dataset or dataarray.
        grid (dict): Grid used for interpolation, dict(lat=[...], lon=[...]).
    """
    lat_min = ds['lat'].min().data if lat_range is None else lat_range[0] 
    lat_max = ds['lat'].max().data if lat_range is None else lat_range[1] 
    lon_min = ds['lon'].min().data if lon_range is None else lon_range[0] 
    lon_max = ds['lon'].max().data if lon_range is None else lon_range[1] 
    init_lat = np.arange(
        lat_min, (lat_max + step_lat), step_lat
    )
    init_lon = np.arange(
        lon_min, lon_max, step_lon
    )
    grid = {'lat': init_lat, 'lon': init_lon}
    # Interpolate
    da = ds.interp(grid, method='nearest')

    return da, grid


def interp_points(i, da, points_origin, points_grid):
    """Interpolation of dataarray to a new set of points.

    Args:
        i (int): Index of time
        da (xr.Dataarray): Dataarray
        points_origin (np.ndarray): Array of origin locations.
        points_grid (np.ndarray): Array of locations to interpolate on.

    Returns:
        i (int): Index of time
        values_grid_flat (np.ndarray): Values on new points.
    """
    values_origin = da[i].data.flatten()
    values_grid_flat = interpolate.griddata(
        points_origin, values_origin, xi=points_grid, method='nearest'
    )
    return i, values_grid_flat


def interp_points2mercato(da, grid, n_cpus=1):
    """Interpolate Dataarray with non-rectangular grid to mercato grid.

    Args:
        da (xr.Dataarray): Dataarray with non-rectangular grid. 
        grid (dict): Grid to interpolate on dict(lat=[...], lon=[...]).

    Returns:
        da_grid (xr.Dataarray): Dataarray interpolated on mercato grid. 
    """
    print(f"Interpolate data with non-rectangular grid to mercato grid. n_cpus={n_cpus}.",
          flush=True)
    # Create array of points from mercato grid
    xx, yy = np.meshgrid(grid['lon'], grid['lat'])
    points_grid = np.array([xx.flatten(), yy.flatten()]).T
    points_origin = np.array(
        [da['nav_lon'].data.flatten(), da['nav_lat'].data.flatten()]).T

    # Interpolation at each time step in parallel
    n_processes = len(da['time_counter'])
    results = Parallel(n_jobs=n_cpus)(
        delayed(interp_points)(i, da, points_origin, points_grid)
        for i in tqdm(range(n_processes))
    )
    # Read results
    ids = []
    values_grid_flat = []
    for r in results:
        i, data = r
        ids.append(i)
        values_grid_flat.append(data)
    ids = np.array(ids)
    values_grid_flat = np.array(values_grid_flat)

    # Store to new dataarray
    values_grid = np.reshape(
        values_grid_flat,
        newshape=(len(values_grid_flat), len(grid['lat']), len(grid['lon']))
    )
    times = da['time_counter'].data[ids]
    da_grid = xr.DataArray(
        data=values_grid,
        dims=['time', 'lat', 'lon'],
        coords=dict(time=times, lat=grid['lat'], lon=grid['lon']),
        name=da.name
    )
    return  da_grid


def cut_map(ds, lon_range=None, lat_range=None, shortest=True):
    """Cut an area in the map. Use always smallest range as default.
    It lon ranges accounts for regions (eg. Pacific) that are around the -180/180 region.

    Args:
    ----------
    lon_range: list [min, max]
        range of longitudes
    lat_range: list [min, max]
        range of latitudes
    shortest: boolean
        use shortest range in longitude (eg. -170, 170 range contains all points from
        170-180, -180- -170, not all between -170 and 170). Default is True.
    Return:
    -------
    ds_area: xr.dataset
        Dataset cut to range
    """
    if lon_range is not None:
        if (max(lon_range) - min(lon_range) <= 180) or shortest is False:
            ds = ds.sel(
                lon=slice(np.min(lon_range), np.max(lon_range)),
                lat=slice(np.min(lat_range), np.max(lat_range))
            )
        else:
            # To account for areas that lay at the border of -180 to 180
            ds = ds.sel(
                lon=ds.lon[(ds.lon < min(lon_range)) |
                           (ds.lon > max(lon_range))],
                lat=slice(np.min(lat_range), np.max(lat_range))
            )
    if lat_range is not None:
        ds = ds.sel(
            lat=slice(np.min(lat_range), np.max(lat_range))
        )

    return ds


def select_months(ds, months=[12, 1, 2]):
    """Select only some months in the data.

    Args:
        ds ([xr.DataSet, xr.DataArray]): Dataset of dataarray
        months (list, optional): Index of months to select. 
                                Defaults to [12,1,2]=DJF.
    """
    ds_months = ds.sel(time=np.in1d(ds['time.month'], months))

    return ds_months


def select_time_snippets(ds, time_snippets):
    """Cut time snippets from dataset and concatenate them.
ra
    Parameters:
    -----------
    time_snippets: np.datetime64  (n,2)
        Array of n time snippets with dimension (n,2).

    Returns:
    --------
    xr.Dataset with concatenate times
    """
    ds_lst = []
    for time_range in time_snippets:
        ds_lst.append(ds.sel(time=slice(time_range[0], time_range[1])))

    ds_snip = xr.concat(ds_lst, dim='time')

    return ds_snip


def average_time_periods(ds, time_snippets):
    """Select time snippets from dataset and average them.

    Parameters:
    -----------
    time_snippets: np.datetime64  (n,2)
        Array of n time snippets with dimension (n,2).

    Returns:
    --------
    xr.Dataset with averaged times
    """
    ds_lst = []
    for time_range in time_snippets:
        temp_mean = ds.sel(time=slice(time_range[0], time_range[1])).mean('time')
        temp_mean['time'] = time_range[0] +  0.5 * (time_range[1] - time_range[0])
        ds_lst.append(temp_mean)

    ds_snip = xr.concat(ds_lst, dim='time')

    return ds_snip


def get_mean_time_series(da, lon_range, lat_range, time_roll=0):
    """Get mean time series of selected area.

    Parameters:
    -----------
    da: xr.DataArray
        Data
    lon_range: list
        [min, max] of longitudinal range
    lat_range: list
        [min, max] of latiduninal range
    """
    da_area = cut_map(da, lon_range, lat_range)
    ts_mean = da_area.mean(dim=('lon', 'lat'), skipna=True)
    ts_std = da_area.std(dim=('lon', 'lat'), skipna=True)
    if time_roll > 0:
        ts_mean = ts_mean.rolling(time=time_roll, center=True).mean()
        ts_std = ts_std.rolling(time=time_roll, center=True).mean()

    return ts_mean, ts_std




def normalize(da, method='zscore'):
    """Normalize dataarray by a given method.

    Args:
        da ([type]): [description]
        method (str, optional): Normalization method. 'minmax' corresponds to 0-1,
            and 'zscore' standardizes the data. Defaults to 'zscore'.

    Returns:
        [type]: [description]
    """
    print(f'Normalize data by {method}!')
    flatten = da.stack(z=da.dims)
    if method == 'minmax':
        min = flatten.min(skipna=True)
        max = flatten.max(skipna=True)
        norm_data = (
            (flatten - min) / (max - min)
        )
        attr = dict(norm=method, min=min.data, max=max.data)
    elif method == 'zscore':
        mean = flatten.mean(skipna=True)
        std = flatten.std(skipna=True)
        norm_data = (
            (flatten - mean) / std
        )
        attr = dict(norm=method, mean=mean.data, std=std.data)
    else:
        print(f'Your selected normalization method "{method}" does not exist.')

    norm_data = norm_data.unstack('z')
    for key, val in attr.items():
        norm_data.attrs[key] = val

    return norm_data


def unnormalize(dmap, attr):
    """Unnormalize data.

    Args:
        dmap (xr.Dataarray): Datamap.
        attr (dict): Dictionary containing normalization information
            attr = {'norm': 'minmax' ,'min': , 'max': } 
            or attr = {'norm': 'zscore' ,'mean': , 'std': } 

    Returns:
        rec_map (xr.Dataarray): Unnormalized map.
    """
    if attr['norm'] == 'minmax':
        rec_map = dmap * (attr['max'] - attr['min']) + attr['min']
    elif attr['norm'] == 'zscore':
        rec_map = (dmap * attr['std'] + attr['mean'])
    else:
        print(
            f'Your selected normalization method {attr["norm"]} does not exist.')

    return rec_map


def rotate_matrix(M, Theta):
    """Rotate 2d matrix by angle theta.

    Args:
        M (np.ndarray): (2,2) 2d matrix
        Theta (float): Angle in rad.

    Returns:
        (np.ndarray) (2,2) Rotated matrix.
    """
    R = np.array(
        [[np.cos(Theta), -np.sin(Theta)], [np.sin(Theta), np.cos(Theta)]]
    )
    return R @ M @ R.T
