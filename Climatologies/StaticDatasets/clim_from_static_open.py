# It computes the open sea seasonal and annual climatology of MedBGCins variables
# Open sea is defined by the 200 m depth and profiles shallower than 50 m are neglected
# Depth layers are the same considered for Emodnet_int2018 dataset
# Also statistics are the same, except for annual mean of N1p, N3n and N5s
# for which interpolation by a logistic curve has been implemented.
# A NetCDF with the same structure of the previous clim is provided for each variable,
# but with added nprofiles and explicit global attributes

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons import season,timerequestors
from bitsea.static.Nutrients_reader import NutrientsReader
from bitsea.static.Carbon_reader import CarbonReader
from bitsea.instruments.var_conversions import NUTRVARS
import numpy as np
from bitsea.commons.layer import Layer
import scipy.optimize as opt
from config import LayerList, REQUESTORS_LIST, basV2
from bitsea.commons.mask import Mask

# EDIT part
modelvarname = "N1p"
# end EDIT part

def corestatistics(pres, values, times):
    STATS = np.zeros((10,),np.float32)
    STATS[0]   = values.mean()
    STATS[1]   = values.std()
    STATS[2:7] = np.percentile(values, [1, 25, 50, 75, 99] )
    STATS[7]   = pres.mean()
    STATS[8]   = times.mean()
    STATS[9]   = len(values)
    return STATS

def mylogistic4(xdata,*x4):
    a = x4[0]
    b = x4[1]
    c = x4[2]
    d = x4[3]
    F = d + (a - d) / (1+ pow((xdata / c),b) )
    return F

def DatasetInfo(modelvarname):
    '''
    1) Helps the user to navigate inside static dataset,
    because some variables is in Nutrients one, some other
    are in the Carbon one.
    2) Defines the official name of variables to be used in deliverables
       statistics.

    Argument:
    * modelvarname * string, like 'N1p'
    Returns:
    * var     * variable name to access dataset
    * dataset * NutrientsReader or Carbonreader object
    '''
    if modelvarname in ['N1p','N3n','O2o','N4n','N5s','P_l','N2n']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    if modelvarname in ['O3h', 'Ac', 'ALK'] :
        var = 'ALK'
        dataset = C
    if modelvarname in ['O3c', 'DIC'] :
        var ='DIC_merged'
        dataset = C
    if modelvarname in ['pH', 'PH'] :
        var='pH_ins_merged'
        dataset = C
    if modelvarname == 'pCO2' :
        var='pCO2_rec'
        dataset = C
    if modelvarname == 'salinity' :
        var='salinity'
        dataset = C
    if modelvarname == 'temp' :
        var='temp'
        dataset = C
    if modelvarname not in ['N1p','N3n','O2o','N4n','N5s','O3h', 'Ac','ALK','O3c', 'DIC', 'pH',',PH', 'pCO2','P_l','N2n', 'salinity','temp']:
        raise ValueError("variable not in static dataset ")
    return var, dataset

N=NutrientsReader()
C=CarbonReader()

var, Dataset = DatasetInfo(modelvarname)

if modelvarname=='N1p':
    useLogistic = True
    startpt=np.asarray([0.1, 0.1, 1000, 0.4],dtype=np.float64)
elif ((modelvarname=='N3n') | (modelvarname=='N5s')): 
    useLogistic = True
    startpt = np.asarray([0.1, 0.1, 500, 4],dtype=np.float64)
else: # default values
    startpt=np.asarray([0.1, 0.1, 1000, 0.4],dtype=np.float64) 
    useLogistic= False

# References for the climatology
TI = TimeInterval('19950101','20240101',"%Y%m%d")
meshmaskfile ='/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc'
TheMask= Mask(meshmaskfile)
mask200 = TheMask.mask_at_level(200.0)
obsdepthlim=50.0
#

nSub    = len(basV2.P.basin_list)
nSeas   = len(REQUESTORS_LIST)
nLayers = len(LayerList)
nStat   = 10

SubDescr    = "" 
DepthDescr  = ""
StatDescr   = "Mean, Std, p01, p25, p50, p75, p99, depth_mean, month_mean, nvalues"
SeasDescr   = ""

for i in basV2.P.basin_list : SubDescr   +=str(i.name) + ", "
for iS in REQUESTORS_LIST : SeasDescr    +=str(iS.string) + ", "

depthLayer = np.zeros(nLayers, np.float32)*np.nan
for ilayer, layer in enumerate(LayerList):
    depthLayer[ilayer] = 0.5*(layer.top+layer.bottom)
    DepthDescr +=str(layer.bottom) + ", "

CLIM = np.zeros((nSeas,nSub,nLayers,nStat), np.float32)*np.nan
NPROFS  = np.zeros((nSeas,nSub), np.int32)

# Checks on the possible logistic computation
if useLogistic==True:
    if startpt.shape!=(4,):
        raise ValueError("the first guess for parameters is not an array of 4 values, as requested ")
    if startpt.dtype!=np.float64:
        raise ValueError("the first guess for parameters is not an array of float64, as requested ")
    if modelvarname not in ['N1p','N3n','N5s']:
        raise ValueError("variable not properly interpolable by a logistic curve (as instead N1p, N3n, N5s)")

for iSeas, seasreq in enumerate(REQUESTORS_LIST):
    print(seasreq.string)
    if seasreq.string=='annual':
        T_int = TI
    else:
        S=season.season()
        T_int = timerequestors.Clim_season(iSeas,S)
    for isub, sub in enumerate(basV2.P):
        print(sub)
        ProfilelistAll = Dataset.Selector(var, T_int, sub) 
        Profilelist=[p for p in ProfilelistAll if TI.contains(p.time)]
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        Months=np.zeros((0,),np.float32)
        nProfiles=0 
        for p in Profilelist:
            lonp = p.lon
            latp = p.lat
            ilon,ilat = TheMask.convert_lon_lat_to_indices(lonp,latp)
            if mask200[ilat,ilon]:
                pres, profile, _ = p.read(var)
                if max(pres)>=obsdepthlim: # to delete observations that are only in the surface layers (not real profiles)
                    Pres   = np.concatenate((Pres,pres))
                    Values = np.concatenate((Values,profile))
                    Months = np.concatenate((Months, np.ones_like(pres)*p.time.month))
                    nProfiles=nProfiles+1
        NPROFS[iSeas,isub] = nProfiles            
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<layer.bottom)
            if (ii.sum()> 1 ) :
                local_pres = Pres[ii]
                local_profile = Values[ii]
                local_month = Months[ii]
                CLIM[iSeas,isub,ilayer,:] = corestatistics(local_pres, local_profile, local_month)
        # if useLogistic=True and case of annual clim, the CLIM mean is overwritten by using results from logistic curve down to 1000 m
        if (useLogistic and seasreq.string=='annual'):
            #print("I try to compute the logistic in ",sub)
            Pres_fit = Pres[Pres<=1000]
            Values_fit = Values[Pres<=1000]
            if nProfiles>0:
                try:
                    popt, pcov = opt.curve_fit(mylogistic4,Pres_fit,Values_fit,p0=startpt,method="lm")
                    #print("popt: ",popt)
                    clim = mylogistic4(depthLayer, *popt)
                except RuntimeError:
                    print("Error - curve_fit failed")
                    clim=np.nan*np.ones_like(depthLayer)
            else:
                clim=np.nan*np.ones_like(depthLayer)
            CLIM[iSeas,isub,:,0] = clim

outfile = "/g100_scratch/userexternal/vdibiagi/clim_MedBGCins/NEW/clim." + modelvarname + ".nc"

import scipy.io as NC
ncOUT = NC.netcdf_file(outfile,"w")
ncOUT.createDimension("seas", nSeas)
ncOUT.createDimension("nSub", nSub)
ncOUT.createDimension("layers", nLayers)
ncOUT.createDimension("Statistics", nStat)
ncvar = ncOUT.createVariable("CLIMATOLOGY", "f", ('seas','nSub','layers','Statistics'))
ncvar[:] = CLIM
ncvar = ncOUT.createVariable("NPROFILES", "f", ('seas','nSub'))
ncvar[:] = NPROFS
setattr(ncOUT,"seas__list"  ,     SeasDescr[:-2])
setattr(ncOUT,"sub___list"  ,     SubDescr[:-2])    
setattr(ncOUT,"PresDOWN_list"  ,  DepthDescr[:-2])
setattr(ncOUT,"stat__list"  ,     StatDescr    )

ncOUT.close()


import pylab as pl
PresCEN = [(l.bottom+l.top)/2  for l in LayerList]
outdir = "/g100_scratch/userexternal/vdibiagi/clim_MedBGCins/NEW/PLOTS/"

for it, seasreq in enumerate(REQUESTORS_LIST):
    for isub, sub in enumerate(basV2.P):
        
        Data_for_plot = CLIM[it,isub,:,:] #nStat x nLayers
        if np.all(np.isnan(Data_for_plot)) :
            print("No data for ", sub.name, seasreq.string)
            continue

        fig = pl.figure(figsize=(30,10))
        ax1=pl.subplot(1,5,1)
        p1, =pl.plot(Data_for_plot[:,0],PresCEN,'k-', label='mean')
        p2, =pl.plot(Data_for_plot[:,2],PresCEN,'b--',label='p01')
        p3, =pl.plot(Data_for_plot[:,3],PresCEN,'b-', label='p25')
        p4, =pl.plot(Data_for_plot[:,4],PresCEN,'k--',label='med')
        p5, =pl.plot(Data_for_plot[:,5],PresCEN,'b-' ,label='p75')
        p6, =pl.plot(Data_for_plot[:,6],PresCEN,'b--',label='p99')
        
        
        pl.legend(handles=[p1,p2,p3,p4,p5,p6],loc=3)
        pl.title(modelvarname +  ' ' + seasreq.string)

        ax2=pl.subplot(1,5,2)
        p1, =pl.plot(Data_for_plot[:,1],PresCEN,'k-',label='std')
        p2, =pl.plot(Data_for_plot[:,5]-Data_for_plot[:,3],PresCEN,'k--',label='irq')
        pl.legend(handles=[p1,p2])

        ax3=pl.subplot(1,5,3)
        p1, =pl.plot(Data_for_plot[:,9],PresCEN,'k*-',label='#obs')
        pl.legend(handles=[p1])

        ax4=pl.subplot(1,5,4)
        p1, =pl.plot(np.ones(nLayers),Data_for_plot[:,7],'k*',label='Zobs')
        p2, =pl.plot(np.ones(nLayers),PresCEN,'r+',label='Zref')
        pl.legend(handles=[p1,p2])

        ax5=pl.subplot(1,5,5)
        p1, =pl.plot(Data_for_plot[:,8],Data_for_plot[:,7],'k*',label='month')
        pl.legend(handles=[p1])

        for ax in [ax1, ax2, ax3, ax4, ax5] : ax.invert_yaxis()
        
        outfilename = outdir + var + '_' + sub.name + '_' + seasreq.string + '.png'
        fig.savefig(outfilename, dpi=150, facecolor='w', edgecolor='w',orientation='landscape')
        pl.close(fig)



