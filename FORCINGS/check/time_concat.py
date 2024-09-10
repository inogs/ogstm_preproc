from bitsea.commons import netcdf4
from bitsea.commons.Timelist import TimeList
import netCDF4
import numpy as np
from bitsea.commons.mask import Mask
import datetime

var="KE_ratio"
outfile="pippo.nc"
INPUTDIR="/g100_scratch/userexternal/gbolzon0/BI-HOURLY/2H/metrics_2d/"
TheMask = Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask_CMCC.nc', loadtmask=False)
jpk,jpj, jpi  = TheMask.shape

TL=TimeList.fromfilenames(None, INPUTDIR, "metrics.2019030*", prefix="metrics.")
Dref = datetime.datetime(1970,1,1,0,0,0)

nFrames = TL.nTimes

A = np.zeros((nFrames,jpj, jpi))
TIME = np.zeros((nFrames,),np.double)

for iframe, filename in enumerate(TL.filelist):
    A[iframe,:] = netcdf4.readfile(filename, var)
    D = TL.Timelist[iframe] - Dref
    TIME[iframe] = D.days*3600*24 + D.seconds

ncOUT = netCDF4.Dataset(outfile,'w')
ncOUT.createDimension('time',nFrames)
ncOUT.createDimension('longitude',jpi)
ncOUT.createDimension('latitude',jpj)

ncvar = ncOUT.createVariable('time','d',('time',))
setattr(ncvar,'units',       'seconds since 1970-01-01 00:00:00')
setattr(ncvar,'long_name'    ,'time')
setattr(ncvar,'standard_name','time')
setattr(ncvar,'axis'         ,'T')
setattr(ncvar,'calendar'     ,'standard')
ncvar[:] = TIME

ncvar[:]

ncvar=ncOUT.createVariable(var,'f',('time','latitude','longitude'), zlib=False, fill_value=1.0e+20)
ncvar[:] = A
ncOUT.close()


import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
i,j = 865,164
fig,ax=pl.subplots()
ax.plot(A[:,j,i])
fig.show()