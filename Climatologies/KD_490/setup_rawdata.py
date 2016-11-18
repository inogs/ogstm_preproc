import netCDF4
import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import Sat.SatManager as Sat
import dom_dec

maskfile="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/KD490_1km_meshmask.nc"
INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"
TI = TimeInterval("19500101","20500101","%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')
nFrames = TL.nTimes

D=netCDF4.Dataset(maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 2
nproc_j = 2

for jj in range(nproc_j):
    for ji in range(nproc_i):
        outfile = "Kd_raw_data_%d_%d" %(jj,ji)
        print outfile
        #if ((jj==0) & (ji==0)) : continue
        
        I_start,I_end, J_start, J_end = dom_dec.dom_dec(ji, jj, jpiglo, jpjglo, nproc_i, nproc_j)
        tmask = tmask_glo[J_start:J_end, I_start:I_end]
        print I_start,I_end, J_start, J_end

        Nwaterpoints = tmask.sum()
        print "Nwaterpoints", Nwaterpoints
        
        RAW_DATA = np.zeros((Nwaterpoints,nFrames),np.float32)
        for ifile, filename in enumerate(TL.filelist):
            A=Sat.readfromfile(filename,'KD490')
            A=A[J_start:J_end, I_start:I_end]
            RAW_DATA[:,ifile] = A[tmask]
        
        np.save(outfile,RAW_DATA)







