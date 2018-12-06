import numpy as np
from commons.genUserDateList import getTimeList
from commons.Timelist import TimeList
from commons import timerequestors
from commons import netcdf4
from commons.mask import Mask
from analysis_config import maskfile, VARLIST, Seas_obj, mydtype


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
nSeas = Seas_obj.numbers_season
nVar = len(VARLIST)


dard_slice = np.load("dard_slice.npy") #np.zeros((14,3,nFrames),np.float32)


MODEL_AEG = np.load("AEG_integrals_24.npy")


LUDWIG_INPUT_KTy={"N1p": 0.2, "N3n":13, "N5s":13, "O3c":3363.8, "O3h": 360.25}
LUDWIG_INPUT={}
w= 1.0e+12;
t = 1./(365 * 86400)
n = 1./14
p = 1./31
s = 1./28
cn = w*n
cp = w*p
cs = w*s
ca = w  
cc = w    
LUDWIG_INPUT['N1p'] =LUDWIG_INPUT_KTy["N1p"]*cp # mmol/y
LUDWIG_INPUT['N3n'] =LUDWIG_INPUT_KTy["N3n"]*cn # mmol/y
LUDWIG_INPUT['N5s'] =LUDWIG_INPUT_KTy["N5s"]*cs # mmol/y
LUDWIG_INPUT['O3h'] =LUDWIG_INPUT_KTy["O3h"]*ca #mg/y
LUDWIG_INPUT['O3c'] =LUDWIG_INPUT_KTy["O3c"]*cc #mg/y

output_nseasons = 1
BOUNDARY_CONCENTRATION=np.zeros((1,),dtype=mydtype)

for var in VARLIST:
    outlet_flux = 0
    for iSeas in range(nSeas):
        req = timerequestors.Season_req(2017,iSeas,Seas_obj)
        #req=timerequestors.Clim_season(iSeas,Seas_obj)
        ii,w = TL.select(req)
        
        Good = np.zeros((nFrames),np.bool)
        for k in ii : Good[k] = True
        u_mean =dard_slice[:,:,Good].mean(axis=2)

        LOWER_flowrate = (u_mean[lev_change:,:] * area[lev_change:,:]).sum() # m3/s
        lower_conc = MODEL_AEG[var][iSeas,1] #mmol/m3
        outlet_flux += LOWER_flowrate*lower_conc*req.time_interval.length()  #mmol (/y)  
    
    inlet_annual_flux = (LUDWIG_INPUT[var] + outlet_flux )
    
    if output_nseasons == 1:
        u_annual_mean = dard_slice.mean(axis=2)
        UPPER_flowrate = (-u_annual_mean[:lev_change,:] * area[:lev_change,:]).sum()
        boundary_conc = inlet_annual_flux/(UPPER_flowrate*86400*365)
        BOUNDARY_CONCENTRATION[var][0]  = boundary_conc


np.save("dard_boundary_values.npy",BOUNDARY_CONCENTRATION)