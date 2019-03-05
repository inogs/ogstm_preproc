from commons.mask import Mask
from commons.dataextractor import DataExtractor
import numpy as np
from commons.Timelist import TimeInterval, TimeList
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1


TheMask=Mask("/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc")
jpk, jpj, jpi = TheMask.shape
E1T=np.zeros((jpk,jpj,jpi),np.float32)
E2T=np.zeros((jpk,jpj,jpi),np.float32)
for k in range(jpk):
    E1T[k,:,:] = TheMask.e1t
for k in range(jpk):
    E2T[k,:,:] = TheMask.e2t


INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_07/wrkdir/MODEL/FORCINGS/"
OUTPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/monthly_DeltaT/"

TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", filtervar="U", prefix="U",hour=0)
eps=1.e-08
Cmax=1.0

MONTHLY_REQ=TL.getMonthlist()
for req in MONTHLY_REQ[rank::nranks]:
    outfile=OUTPUTDIR + req.string + ".txt"
    ii,w = TL.select(req)
    nFrames=len(ii)
    DELTAT=np.zeros((nFrames,5),np.float32)
    
    for iframe, k in enumerate(ii):
        timestr = TL.Timelist[k].strftime("%Y%m%d-%H:%M:%S")
        print timestr
        filenameU=INPUTDIR + "U" + timestr + ".nc"
        filenameV=INPUTDIR + "V" + timestr + ".nc"
        filenameW=INPUTDIR + "W" + timestr + ".nc"

        U=np.abs(DataExtractor(TheMask,filenameU,"vozocrtx").values)
        V=np.abs(DataExtractor(TheMask,filenameV,"vomecrty").values)
        W=np.abs(DataExtractor(TheMask,filenameW,"vovecrtz").values)
        
        U[U==0]=eps
        V[V==0]=eps
        W[W==0]=eps
        Fact = U/E1T + V/E2T + W/TheMask.e3t
        deltat = Cmax/Fact#[:,70:170,200:382]
        low=deltat <450
        K,J,I = np.nonzero(deltat==deltat.min())
        DELTAT[iframe,0]=deltat.min()
        DELTAT[iframe,1]=low.sum()
        DELTAT[iframe,2]=K[0]
        DELTAT[iframe,3]=J[0]
        DELTAT[iframe,4]=I[0]
    np.savetxt(outfile, DELTAT,fmt="%10.3f %d %d %d %d")
    


