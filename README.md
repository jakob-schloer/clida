# clida

This is a collection of scripts for downloading weather and climate data and preprocessing them accordingly.

## Set-up

### Required packages

The following packages are required for running the scripts in this repo:
- xarray 
- netcdf4
- pandas
- numpy
- scipy
- gcsfs
- zarr
- cdsapi
- cdo
- cftime
- re
- tqdm
- joblib

### Download ERA5

ERA5 data are downloaded via the Climate Data Store (CDS) Application Program Interface (API). Follow the steps on their [how-to-api](https://cds.climate.copernicus.eu/api-how-to).


## How-to use the scripts 

1. Copy the `dwnld_config_template.py` to `dwnld_config.py`
2. Change the `"raw_data_dir"` in the `dwnld_config.py` to your preferred location
3. Select only the part in the config which you need and adapt it to your needs
4. Run your `download_<product>.py` with python. Keep in mind that the `dwnld_config.py` is in the same directory

To run scripts on the slurm queueing system, I've provided a sample script in `example_scripts/submit_to_slurm.sbatch` which can be started using `sbatch`. 

## TODO
- WeatherBench
- CPCP, TRMM

## Contributing 
- all final output files should be saved to `climate_data/[model or data product name]`
- all final output files should be in NetCDF format
- all final output files should be in a one-file-per-variable format, with the variable in the file called the same as the `[variable shorthand]` in the filename (below)

Furthermore, please document your code as well as you can. Additionally, when you can, add links to documentation on how to access the datasets, or any additional tools that are needed (for example, ECMWF's API). 




