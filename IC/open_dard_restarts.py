from commons import netcdf4
from commons.mask import Mask
from IC import RSTwriter
import os,glob
import numpy as np

INPUTDIR="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/IC/RST_CLOSE_DARD/"
OUTPUTDIR="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/IC/RST_OPEN__DARD/"
TheMask=Mask("/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/meshmask.nc")

filelist=glob.glob(INPUTDIR+"*nc")

I=[841,841,841]
J=[235,236,237]
nOpen=len(I)


for inputfile in filelist:

    basename=os.path.basename(inputfile)
    print basename
    var=basename[22:25]
    outfile=OUTPUTDIR + basename
    continue
    R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]
    for k in range(nOpen):
        i=I[k]
        j=J[k]
        R[:,j,i]=R[:,j,i-1]
    RSTwriter(outfile, var, R, TheMask)
    if var == "R7c":
        outfile=outfile.replace("R7c","R3c")
        RSTwriter(outfile, "R3c", R, TheMask)



# O5c is defined on the basis of 3 values, it does not exist in BFMv2
jpk, jpj, jpi = TheMask.shape
PIC = np.ones((jpk,jpj,jpi),np.double)*1.e-11

jk = TheMask.getDepthIndex(300)
PIC[:jk,:,:]=1.0
jk = TheMask.getDepthIndex(150)
PIC[:jk,:,:]=2.0
outfile = outfile.replace(var,"O5c")

RSTwriter(outfile, 'O5c', PIC, TheMask)
