
from commons.mask import Mask
from bclib.bounmask import bounmask
from bclib import atmosphere
from bclib.gib import gib
from bclib.co2 import co2atm

import config as conf
TheMask = Mask(conf.file_mask)

CO2 =co2atm(conf)
CO2.generate(TheMask)
ATM=atmosphere.atmosphere(TheMask,conf)
ATM.write_netcdf(TheMask, conf.dir_out)


BOUN = bounmask(conf)
BOUN.generate(TheMask)
BOUN.write_netcdf(TheMask)

G = gib(conf,TheMask)
G.generate(TheMask, BOUN)





import sys
sys.exit()
## RIVER, to be tested
index = BOUN.load('index')

from bclib.river import river_data
R = river_data(conf)
R.modularize(conf)
georef = R.gen_map_indexes(TheMask)
idxt, positions = R.gen_boun_indexes(index)

R.generate_monthly_files(conf, TheMask,idxt, positions)

climatological=False
if climatological:
    for sheet in conf.river_data_sheet:
        R.river_data[sheet]['yyyy'] = (R.river_data[sheet]['2000'] + R.river_data[sheet]['2001'] + R.river_data[sheet]['2002'])/3
    R.generate_climatological_monthly_files(conf, TheMask)
