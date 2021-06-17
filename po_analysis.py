import numpy as np
from commons import netcdf4
from commons.Timelist import TimeList

INPUTDIR="/gpfs/scratch/userexternal/gcoidess/EAS6_2019_PO/"

riverinput="/gpfs/scratch/userexternal/gbolzon0/PO/runoff_1d_nomask_y2019.nc"

S = netcdf4.readfile(riverinput, 's_river')
R = netcdf4.readfile(riverinput, 'sorunoff')


ii = S[0,:,:]==18
J,I = np.nonzero(ii)

nPoints=ii.sum()

TL = TimeList.fromfilenames(None, INPUTDIR, "*T.nc", prefix="assw_drpo_2019_1d_", dateformat="%Y%m%d")

nFrames = TL.nTimes
PO_RUNOFF=np.zeros((nFrames,nPoints), np.float32)
INPUT_RUNOFF=np.zeros((nFrames,nPoints), np.float32)


for iFrame, filename in enumerate(TL.filelist):
    print iFrame
    nemo_Runoff=netcdf4.readfile(filename, "sorunoff")[0,:]
    PO_RUNOFF[iFrame,:] = nemo_Runoff[ii]
    INPUT_RUNOFF[iFrame,:] = R[iFrame,ii]
    
np.save('nemo_runoff.npy',PO_RUNOFF)
np.save('input_runoff.npy',INPUT_RUNOFF)


import pylab as pl


for k in range(10):
    fig,ax = pl.subplots()
    ax.plot(TL.Timelist,PO_RUNOFF[:,k],'b-',label='ref')
    ax.plot(TL.Timelist,INPUT_RUNOFF[:,k],'r.',label='model')
    ax.legend()
    outfile="P%d.png" %(k)
    fig.savefig(outfile)
    pl.close(fig)

    