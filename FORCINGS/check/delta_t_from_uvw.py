from commons.mask import Mask
from commons.dataextractor import DataExtractor
import numpy as np
TheMask=Mask("/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc")

INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_07/wrkdir/MODEL/FORCINGS/"
timestr="20180415-00:00:00"
filenameU=INPUTDIR + "U" + timestr + ".nc"
filenameV=INPUTDIR + "V" + timestr + ".nc"
filenameW=INPUTDIR + "W" + timestr + ".nc"

U=DataExtractor(TheMask,filenameU,"vozocrtx").values
V=DataExtractor(TheMask,filenameV,"vomecrty").values
W=DataExtractor(TheMask,filenameW,"vovecrtz").values

eps=1.e-08
ii=W==0
W[ii] = eps
W[~TheMask.mask]=eps

ii=U==0
U[ii] = eps
ii=V==0
V[ii] = eps

jpk, jpj, jpi = TheMask.shape
E1T=np.zeros((jpk,jpj,jpi),np.float32)
E2T=np.zeros((jpk,jpj,jpi),np.float32)
for k in range(jpk):
    E1T[k,:,:] = TheMask.e1t
for k in range(jpk):
    E2T[k,:,:] = TheMask.e2t

deltat_w=TheMask.e3t/np.abs(W)
deltat_w[~TheMask.mask]=1.e+20

deltat_u =E1T/abs(U)
deltat_v =E2T/abs(V)


print deltat_u.min(), deltat_v,min(), deltat_w.min()




