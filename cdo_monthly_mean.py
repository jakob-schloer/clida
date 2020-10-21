
import numpy as np
import os, sys
import scipy.stats as stats
from cdo import  *

def monthly_mean(infile, outfile):
    """Monthly mean."""
    if os.path.exists(outfile):
        print("File already exists!")
    else:
        cdo.monmean(options='-f nc', output=outfile, input=infile)
        print('Done! Created {} file.'.format(outfile))
    return


if __name__ == "__main__":
    # Create cdo object
    tempdir='./tmp_cdo/'
    cdo = Cdo()

    dirname = "/home/jakob/climate_data/local/era5/"
    
    # First compute monthly mean every year
    years = np.arange(2008, 2009, 1)
    mon_infiles = []
    for y in years:
        infile = dirname + "era5_sst_{}_merged.nc".format(y)
        outfile = dirname + "era5_sst_{}_mon_mean.nc".format(y)
        monthly_mean(infile, outfile)
        mon_infiles.append(outfile)

    # Mege monthly data 
#    cdo.mergetime(options='-f nc', output=dirname+"era5_sst_mon_mean.nc", input=mon_infiles)
    
    

    
    
