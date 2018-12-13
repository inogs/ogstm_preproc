from commons import netcdf4
from commons.Timelist import TimeList,TimeInterval
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import numpy as np
from commons import timerequestors 

#TheMask=Mask("/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_02/wrkdir/MODEL/meshmask.nc")
Mask_0 = TheMask.cut_at_level(0)
area = TheMask.area
jpk, jpj, jpi= TheMask.shape

INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_02/wrkdir/MODEL/BC"
TI = TimeInterval("2017","2018","%Y")
TL =TimeList.fromfilenames(TI, INPUTDIR, "TIN_2*", prefix="TIN_")
nFrames = TL.nTimes

bool_dtype=[(sub.name,np.bool)    for sub in OGS.P]
sub__dtype=[(sub.name,np.float32) for sub in OGS.P]

SUBM = np.zeros((jpj, jpi), dtype=bool_dtype)
for isub, sub in enumerate(OGS.Pred):
    S=SubMask(sub, maskobject=Mask_0)
    SUBM[sub.name]=S.mask
    SUBM['med']=(SUBM['med'] | SUBM[sub.name])

nSub = len(OGS.P.basin_list)

for var in ["N3n","N1p","O2o","N5s","O3h","O3c"]:
    print var
    OUT = np.zeros((nFrames,),dtype=sub__dtype)
    for iFrame, filename in enumerate(TL.filelist):
        B=netcdf4.readfile(filename, "riv_" + var)
        good = B > -1
        for isub, sub in enumerate(OGS.P):
            localmask=SUBM[sub.name]
            V = B*area
            river_on_sub = V[localmask & good].sum()
            OUT[sub.name][iFrame] = river_on_sub #mmol/s or mg/s
            #if river_on_sub > 0:
            #    print sub.name, river_on_sub
        req = timerequestors.Monthly_req(2017,iFrame+1)
        req.time_interval.length()
    import sys
    sys.exit()

for month in range(1,13):
    
    



