
"""Merge .nc files."""

import numpy as np
import os, sys
import scipy.stats as stats
from cdo import Cdo 


if __name__ == "__main__":
    # Create cdo object
    tempdir='./tmp_cdo/'
    cdo = Cdo()

    year_start = int(sys.argv[1])
    year_end = int(sys.argv[2])

    dirname = "/home/jakob/climate_data/local/era5/"
    years = np.arange(year_start, year_end, 1)
    for y in years:
        months = np.arange(1, 13, 1)
        infiles = []
        for m in months:
            subdirname = dirname + "{}/{:02d}".format(y, m)
            for file in os.listdir(subdirname):
                if file.endswith(".nc"):
                    infiles.append(os.path.join(subdirname, file))

        # Merge files
        mergefile = dirname + "era5_sst_{}_merged.nc".format(y)
        cdo.mergetime(options='-f nc', output=mergefile, input=infiles)
        print('Created merge file '+ mergefile)

        # Monthly mean
        meanfile = dirname + "era5_sst_{}_mon_mean.nc".format(y)
        cdo.monmean(options='-f nc', output=meanfile, input=mergefile)
        print('Created monthly mean {} file.'.format(meanfile))

