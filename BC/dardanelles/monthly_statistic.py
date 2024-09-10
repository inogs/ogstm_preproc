import numpy as np
from bitsea.commons import timerequestors
import pylab as pl
from bitsea.commons.mask import Mask
from analysis_config import maskfile
from bitsea.commons import netcdf4
from bitsea.commons.Timelist import TimeList
from bitsea.commons.genUserDateList import getTimeList
I=841
area=np.zeros((14,3),np.float32)


TheMask = Mask(maskfile, dzvarname="e3t_0")
e1t = netcdf4.readfile(maskfile, "e1t")[0,0,:,:]
for k in range(14):
    for jlocal, j in enumerate(range(235,238)):
        area[k,jlocal] = e1t[j,I]*TheMask.e3t[k,j,I]



TL = TimeList(getTimeList("20170102-00:00:00", "20171231-00:00:00", "days=1"))
nFrames = TL.nTimes

lev_change=6 #


dard_slice = np.load("dard_slice.npy") #np.zeros((14,3,nFrames),np.float32)
U_upper=np.zeros((12,), np.float32)
U_lower=np.zeros((12,), np.float32)
for month in range(12):
    req=timerequestors.Monthly_req(2017,month+1)
    ii,w = TL.select(req)
    Good = np.zeros((nFrames),np.bool)
    for k in ii : Good[k] = True
    u_mean =dard_slice[:,:,Good].mean(axis=2)
    UPPER_flowrate = (u_mean[:lev_change,:] * area[:lev_change,:]).sum()
    LOWER_flowrate = (u_mean[lev_change:,:] * area[lev_change:,:]).sum() # m3/s
    U_upper[month] = UPPER_flowrate/area[:lev_change,:].sum()
    U_lower[month] = LOWER_flowrate/area[lev_change:,:].sum()


fig,ax=pl.subplots()
ax.plot(U_upper,'r',label='upper')
ax.plot(U_lower,'b',label='lower')
ax.grid()
ax.legend()
ax.set_xlabel("month")
ax.set_ylabel("m/s")
ax.set_title("U monthly average")
fig.savefig("dardanelles_monthly_average.png")