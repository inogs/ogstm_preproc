from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.commons import netcdf4
import numpy as np


INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS"

TI = TimeInterval("2017","2018","%Y")
TL = TimeList.fromfilenames(TI, INPUTDIR, "*U.nc", prefix="simdd_v5_1d_20170401_", dateformat="%Y%m%d", hour=0)

nFrames = TL.nTimes
dard_slice = np.zeros((14,3,nFrames),np.float32)

for iframe, filename in enumerate(TL.filelist):
    dard_slice[:,:,iframe] = netcdf4.readfile(filename, "vozocrtx")[0,:14,235:238,1062]

np.save('dard_slice',dard_slice)