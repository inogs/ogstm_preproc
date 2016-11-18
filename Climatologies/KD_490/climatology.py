import netCDF4
import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import dom_dec
import time_manager

maskfile="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/KD490_1km_meshmask.nc"
INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"


D=netCDF4.Dataset(maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 2
nproc_j = 2

for jproc in range(nproc_j):
    for iproc in range(nproc_i):
        #if ((jj==0) & (ji==0)) : continue
        inputfile = "Kd_raw_data_%d_%d.npy" %(jproc,iproc)
        outputfile= "Kd_clim_%d_%d" %(jproc,iproc)
        print inputfile
        RAW_DATA = np.load(inputfile)
        
        
        I_start,I_end, J_start, J_end = dom_dec.dom_dec(iproc, jproc, jpiglo, jpjglo, nproc_i, nproc_j)
        tmask = tmask_glo[J_start:J_end, I_start:I_end]
        print I_start,I_end, J_start, J_end

        Nwaterpoints = tmask.sum()
        WaterPoints = np.zeros_like(tmask, dtype=np.int32)
        WaterPoints[tmask] = np.arange(Nwaterpoints)
        print "Nwaterpoints", Nwaterpoints
        jpi = I_end-I_start
        jpj = J_end-J_start
        CLIM = np.zeros((jpj,jpi,366), dtype=[('NUMB',np.int32), ('MEAN',np.float32),('STD',np.float32)])
        for julian in range(366):
            print julian
            II, filelist=time_manager.getfilelist(julian)
            raw_data_julian=RAW_DATA[:,II]
            for ji in range(jpi):
                for jj in range(jpj):
                    if tmask[jj,ji]:
                        jw = WaterPoints[jj,ji];
                        PILAloc = raw_data_julian[jw,:];
                        valid_points =PILAloc>-999
                        if valid_points.sum() < 5: continue                        
                        pilarray= PILAloc[valid_points];

                        m = pilarray.mean()
                        s = pilarray.std()
                        outsiders = pilarray > m+3.0*s
                        pilarray=pilarray[~outsiders]
                        n = len(pilarray)
                        CLIM['NUMB'][jj,ji,julian] = n
                        if n>4:
                            m = pilarray.mean()
                            CLIM['MEAN'][jj,ji,julian] = m
                            CLIM['STD' ][jj,ji,julian] =np.sqrt( ((pilarray-m)**2).sum()/(n-1)  )
                        
        np.save(outputfile,CLIM)
        
