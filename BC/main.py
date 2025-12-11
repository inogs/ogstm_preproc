
from bitsea.commons.mask import Mask
from bclib.bounmask import bounmask
from bclib import atmosphere
from bclib.gib import gib
from bclib.co2 import co2atm
from bclib.river import river
import config as conf
import numpy as np
from dardanelles import dardanelles_generator

TheMask = Mask.from_file(conf.file_mask)

CO2 =co2atm(conf)
CO2.generate(TheMask, experiment="RCP85")
ATM=atmosphere.atmosphere(TheMask,conf)
ATM.write_netcdf(TheMask, conf.dir_out, "Area")


BOUN = bounmask(conf)
BOUN.generate(TheMask)     # can be commented in case of test
BOUN.write_netcdf(TheMask, all_variables=False)

G = gib(conf,TheMask)
G.generate(TheMask, BOUN, all_variables=False)

#R = river(conf) # here excel is read
#R.modularize(conf)
#R.gen_map_indexes(TheMask)
#
#climatological=True
#
#if climatological:
#    YEARS = np.arange(2010,2021)
#    for sheet in conf.river_data_sheet:
#        SUM = np.zeros_like(R.river_data[sheet]['2000'], np.float32)
#        for year in YEARS:
#            year_str = str(year)
#            SUM +=R.river_data[sheet][year_str]
#        R.river_data[sheet]['yyyy'] = SUM/len(YEARS)
#    R.generate_climatological_monthly_files(conf, TheMask)
#else:
#    R.generate_monthly_files(conf, TheMask)
#
#dardanelles_generator.dump_files(TheMask, conf.dir_out)

