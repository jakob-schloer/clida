''' File description

@Author  :   Jakob Schlör 
@Time    :   2022/08/08 17:06:08
@Contact :   jakob.schloer@uni-tuebingen.de
'''
# %%
# Packages 
# ======================================================================================



# Stored paths
links = [
    {'name': 'SODA3.12.2', 'url': 'wget -r -l1 --no-parent --progress=bar -nd -A.nc https://dsrs.atmos.umd.edu/DATA/soda3.12.2/REGRIDED/ocean/'},
    {'name': 'HadISST', 'url': 'wget https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz'},
    {'name': 'COBE2', 'url': 'wget https://downloads.psl.noaa.gov/Datasets/COBE2/sst.mon.mean.nc'},
    {'name': 'ERSSTv5', 'url': 'wget https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc'},
    {'name': 'GODAS', 'url': 'wget https://downloads.psl.noaa.gov/Datasets/godas/pottmp.1980.nc'}, # own script required
    {'name': 'Tropflux', 'url': 'wget https://incois.gov.in/tropflux/DataDownload.jsp'}, # own script required
]


# %%
# Packages 
# ======================================================================================