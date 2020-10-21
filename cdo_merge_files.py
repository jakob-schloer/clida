"""Merge .nc files."""

import numpy as np
import os, sys
import scipy.stats as stats
from cdo import Cdo 



def merge_nc_files(in_directory, years, outfile):
    """Merge nc files stored in a folder dirname which are sorted by
       years/month folder structure.

    Args:   
        passdirname = "/home/jakob/climate_data/local/era5/"
        years       = [2019]
        outfile     = "/home/jakob/climate_data/local/era5/merged_sst_era5.nc"
       
    """

    for y in years:
        months = np.arange(1, 13, 1)
        infiles = []
        for m in months:
            dirname = in_directory + "{}/{:02d}".format(y, m)
            for file in os.listdir(dirname):
                if file.endswith(".nc"):
                    infiles.append(os.path.join(dirname, file))

        outfile = in_directory + "{}/era5_sst_{}_merged.nc".format(y, y)
        cdo.mergetime(options='-f nc', output=outfile, input=infiles)

        print('Created merge file '+ outfile)


# if __name__ == "__main__":
#     # Create cdo object
#     tempdir='./tmp_cdo/'
#     cdo = Cdo()
# 
#     dirname = "/home/jakob/climate_data/local/era5/"
#     years = np.arange(2008, 2019, 1)
#     outfile = dirname 
# 
#     merge_nc_files(dirname, years, outfile)


if __name__ == "__main__":
    # Create cdo object
    tempdir='./tmp_cdo/'
    cdo = Cdo()

    dirname = "/home/jakob/climate_data/local/era5/"
    infiles = [dirname + "era5_sst_2018_mon_mean.nc", dirname + "era5_sst_2019_mon_mean.nc"]

    outfile = dirname + "era5_sst_2018-2019_mon_mean.nc"
    cdo.mergetime(options='-f nc', output=outfile, input=infiles)


