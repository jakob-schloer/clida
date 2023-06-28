# clida

This is a collection of scripts for downloading weather and climate data and preprocessing them accordingly.

## Using these scripts

Specify the output folder directory and parameters in the dwnld_config.py which is imported in each script.

A template can be found in `dwnld_config_template.py`.

The following packages are required for running the scripts in this repo:
- xarray 
- pandas
- numpy
- gcsfs
- zarr
- cdsapi
- cdo


## TODO
- WeatherBench
- CPCP, TRMM

## Contributing 
- all final output files should be saved to `climate_data/[model or data product name]`
- all final output files should be in NetCDF format
- all final output files should be in a one-file-per-variable format, with the variable in the file called the same as the `[variable shorthand]` in the filename (below)

Furthermore, please document your code as well as you can. Additionally, when you can, add links to documentation on how to access the datasets, or any additional tools that are needed (for example, ECMWF's API). 




