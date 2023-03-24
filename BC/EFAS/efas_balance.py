import numpy as np
import river_reader as rr

RIVER_TABLE = rr.RIVERS
nPoints = RIVER_TABLE.size
DISCHARGE = np.load('/g100_scratch/userexternal/gbolzon0/EFAS/preproc/BC/EFAS/efas_discharge.npy')

for ir in range(nPoints):
    daily_discharge_timeseries = DISCHARGE[:,ir]*86400. # m3
    yearly_discharge = daily_discharge_timeseries.sum() # m3
    conc = RIVER_TABLE[ir]['N1p'] #g/m3
    mass = yearly_discharge.sum() * conc #g