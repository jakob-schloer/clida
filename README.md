# clida
This is a collection of scripts for downloading weather and climate data and preprocessing them accordingly.

A lot of these scripts are adopted versions of the climate-download repo of Kevin Schwarzwald and Felix Oertel.
For running scripts on cluster computers all scripts are .py files.`

## Using these scripts

Specify the output folder directory and parameters in the dwnld_config.py which is imported in each script.

A template for CMIP6 data can be found at `dwnld_config_template.py`.

The following packages are required for running the scripts in this repo:
- xarray 
- pandas
- numpy
- gcsfs
- zarr


## Coming soon
- ERA5
- observational sea surface temperature datasets (OISST, ERSST)

## Contributing 
- all final output files should be saved to `climate_data/[model or data product name]`
- all final output files should be in NetCDF format
- all final output files should be in a one-file-per-variable format, with the variable in the file called the same as the `[variable shorthand]` in the filename (below)
- all final output files should follow CMIP5 naming conventions, i.e.: 

`[variable shorthand]_[data frequency]_[model or data product name]_[experiment or timeframe]_[run or other info]_[date in YYYYMMDD-YYYYMMDD format](_[file suffix with geographic information]).nc`

Furthermore, please document your code as well as you can. Additionally, when you can, add links to documentation on how to access the datasets, or any additional tools that are needed (for example, ECMWF's API). 




