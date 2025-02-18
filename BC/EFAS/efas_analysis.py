import numpy as np
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
import river_reader as rr
from bitsea.commons.dataextractor import DataExtractor
import seawater as sw

CMCC_Mask=Mask.from_file('/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask_CMCC.nc')
RIVER_TABLE = rr.RIVERS
nPoints = RIVER_TABLE.size
PRES = np.ones((nPoints,),np.float32)*CMCC_Mask.zlevels[0]
SALI = RIVER_TABLE['SAL']
AREA = np.zeros((nPoints,),np.float32)
for ir in range(nPoints):
    i = RIVER_TABLE[ir]['I'] + 222 -1
    j = RIVER_TABLE[ir]['J'] -1
    AREA[ir] = CMCC_Mask.area[j,i]



INPUTDIR="/g100_scratch/userexternal/gbolzon0/EFAS/EAS6_EFAS/"

TL = TimeList.fromfilenames(None, INPUTDIR, "*T.nc", prefix="simu_EAS6_EFAS_1d_", dateformat="%Y%m%d")
nFrames = TL.nTimes




RUNOFF      = np.zeros((nFrames,nPoints), np.float32)
TEMP        = np.zeros((nFrames,nPoints), np.float32)
DISCHARGE   = np.zeros((nFrames,nPoints), np.float32)

for iFrame, filename in enumerate(TL.filelist):
    print(iFrame)
    Runoff   =  DataExtractor(CMCC_Mask,filename, "sorunoff").values
    Temp      = DataExtractor(CMCC_Mask, filename, 'votemper').values[0,:]
    for ir in range(nPoints):
        i = RIVER_TABLE[ir]['I'] + 222 -1
        j = RIVER_TABLE[ir]['J'] -1
        RUNOFF[iFrame,ir] = Runoff[j,i]
        TEMP[  iFrame,ir] =   Temp[j,i]
        
    T      = sw.temp(SALI,TEMP[iFrame,:],PRES)
    RHO    = sw.dens(SALI,T,PRES)
    
    Discharge = RUNOFF[iFrame,:]*AREA/RHO  # m^3/s
    DISCHARGE[iFrame,:] =  Discharge

np.save('efas_runoff.npy',RUNOFF)
np.save('efas_discharge.npy',DISCHARGE)



