import numpy as np
import river_reader as rr
from bitsea.commons.mask import Mask
from bitsea.basins import V2 as OGS
import os
from bitsea.commons.utils import writetable

# EFAS input
RIVER_TABLE = rr.RIVERS
nPoints = RIVER_TABLE.size
DISCHARGE = np.load('/g100_scratch/userexternal/gbolzon0/EFAS/preproc/BC/EFAS/efas_discharge.npy')


# PERSEUS_INPUT
os.chdir("/g100_scratch/userexternal/vdibiagi/preproc/BC")
import config_efas as conf
from bclib.river import river
conf.file_river = 'Perseus-4.6_39rivers_mesh24.xlsx'
PERSEUS = river(conf)
PERSEUS.modularize(conf)
YEARS = np.arange(2019,2020)
CLIM={}
os.chdir("/g100_scratch/userexternal/vdibiagi/preproc/BC/EFAS")


VARS=['N3n', 'N1p', 'DIC', 'ALK'] # DOC and POC are not in PERSEUS (config with 39 rivers)
nVar=len(VARS)

# EFAS river data are in:
# N3n, N1p, DIC in g/m3, below they will be multiplied for discharge in m3 -> load in g
# ALK Gmol/km3, below it will be multiplied for discharge in m3 -> load in 10^9 mol/10^9 -> load in mol

# Perseus data are
# KM3perYR_NOBLS
# DIP_KTperYR_NOBLS -> 1 g is 10^-9 KT -> conv to compare EFAS data with Perseus data is 10^-9 
# DIN_KTperYR_NOBLS
# DIC_KTperYR_NOBLS
# ALK_GmolperYR_NOBLS -> 1 mol is 10^-9 Gmol -> conv is 10^-9

# from EFAS to Perseus
CONVERSION_DICT={
         'N3n' : 1.0E-09,
         'N1p' : 1.0E-09,
         'DIC' : 1.0E-09,
         'ALK' : 1.0E-09,
         'DOC': 1.0E-09
         }
VARS_SHEETS_DICT={
         'N3n' : 'DIN_KTperYR_NOBLS',
         'N1p' : 'DIP_KTperYR_NOBLS',
         'DIC' : 'DIC_KTperYR_NOBLS',
         'ALK' : 'ALK_GmolperYR_NOBLS'
         }

# We select the different rivers by name, to have a total of 39 rivers 
RIVERS_NAMES = list(dict.fromkeys(RIVER_TABLE['name'])) 
nRiv=len(RIVERS_NAMES)       

# to test the script
# RIVERS_NAMES=['Po']

 
# ANNUAL COMPARISON TABLES
YEARLY_DISCHARGE = np.zeros((nRiv,2),np.float32) 
YEARLY_LOAD = np.zeros((nRiv,nVar,2),np.float32)


for idxr, irn in enumerate(RIVERS_NAMES):

   # EFAS 
   ii=rr.get_indexes_by_river(irn) #all points for the selected river (e.g. Po)
   daily_discharge_timeseries = DISCHARGE[:,ii]*86400. # m3 in each day
   # Here I actually sum in time separately for the points of the same river
   yearly_discharge = daily_discharge_timeseries.sum(axis=0) # m3 in the year
   # Here I sum on river points and use the conv factor to compare data with Perseus (KM3perYR)
   YEARLY_DISCHARGE[idxr,0] = yearly_discharge.sum()*1.0E-09 # total of the river in km3 in the year

   # PERSEUS
   # Find indexes of rivers: from ii to idx1, idx2  as start and end indexes, respectively 
   jj = [i for i, x in enumerate(ii) if x]
   idx1 = jj[0]    
   idx2 = jj[-1]+1  


   os.chdir("/g100_scratch/userexternal/vdibiagi/preproc/BC")

   for idxv, iv in enumerate(VARS):
        
        # EFAS
        VARCONV = CONVERSION_DICT[iv]
        yearly_load_points = RIVER_TABLE[ii][iv]*VARCONV*yearly_discharge # This includes separately all points of the river
        YEARLY_LOAD[idxr, idxv, 0] = yearly_load_points.sum() # sum on the points

        # PERSEUS
        sheet = VARS_SHEETS_DICT[iv]
        SUM = 0.
        SUM_DIS = 0.
        for year in YEARS:
           year_str = str(year)
           annual_contribution=PERSEUS.xls_data[sheet][year_str][idx1:idx2].sum() 
           annual_discharge=PERSEUS.xls_data['KM3perYR_NOBLS'][year_str][idx1:idx2].sum() 
           print(year, sheet, annual_contribution)
           SUM += annual_contribution
           SUM_DIS += annual_discharge
        CLIM[sheet] = SUM/len(YEARS)
        YEARLY_LOAD[idxr, idxv, 1] = SUM/len(YEARS)
        YEARLY_DISCHARGE[idxr,1] = SUM_DIS/len(YEARS) 

####################
OUTPUTDIR='/g100_scratch/userexternal/vdibiagi/articleBC_2023/DATA/'
writetable(OUTPUTDIR + "discharge_2019_total.txt", YEARLY_DISCHARGE, RIVERS_NAMES, ['EFAS','PERSEUS'])

for idxv, iv in enumerate(VARS):
   writetable(OUTPUTDIR + "load_"+ iv + "_2019_total.txt", YEARLY_LOAD[:,idxv,:], RIVERS_NAMES, ['EFAS','PERSEUS'])

 