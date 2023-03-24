import numpy as np
import river_reader as rr
from commons.mask import Mask
from basins import V2 as OGS


RIVER_TABLE = rr.RIVERS
nPoints = RIVER_TABLE.size
DISCHARGE = np.load('/g100_scratch/userexternal/gbolzon0/EFAS/preproc/BC/EFAS/efas_discharge.npy')


for ir in range(nPoints):
    daily_discharge_timeseries = DISCHARGE[:,ir]*86400. # m3
    yearly_discharge = daily_discharge_timeseries.sum() # m3
    conc = RIVER_TABLE[ir]['N1p'] #g/m3
    mass = yearly_discharge * conc #g


TheMask=Mask('/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc',loadtmask=False)
ii = rr.get_indexes_by_river('Po')
ii = rr.get_indexes_by_region(OGS.adr1, TheMask)