# Generates Clim file and images for Static datasets

from static.Nutrients_reader import NutrientsReader
from static.Carbon_reader import CarbonReader

from instruments.var_conversions import NUTRVARS
from config import LayerList, REQUESTORS_LIST, basV2
import numpy as np
N=NutrientsReader()
C=CarbonReader()

def corestatistics(pres, values,times):
    STATS = np.zeros((10,),np.float32)
    STATS[0]   = values.mean()
    STATS[1]   = values.std()
    STATS[2:7] = np.percentile(values, [1, 25,50,75, 99] )
    STATS[7]   = pres.mean()
    STATS[8]   = times.mean()
    STATS[9]   = len(values)
    return STATS

#VARLIST = ['N1p', 'N3n','O2o'] # for nutrient static dataset
modelvarname='N1p'
modelvarname='N3n'
var = NUTRVARS[modelvarname]


nSub    = len(basV2.P.basin_list)
nSeas   = len(REQUESTORS_LIST)
nLayers = len(LayerList)
nStat   = 10

CLIM = np.zeros((nSeas,nSub,nLayers,nStat), np.float32)*np.nan



for iSeas, seasreq in enumerate(REQUESTORS_LIST):
    print seasreq.string
    for isub, sub in enumerate(basV2.P):
        print " --- * ", sub.name
        Profilelist =N.Selector(var, seasreq, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        Months=np.zeros((0,),np.float32)
        for p in Profilelist: 
            pres, profile, Qc = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
            Months = np.concatenate((Months, np.ones_like(pres)*p.time.month))

        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<=layer.bottom)
            if (ii.sum()> 1 ) :
                local_pres = Pres[ii]
                local_profile = Values[ii]
                local_month = Months[ii]
                CLIM[iSeas,isub,ilayer,:] =corestatistics(local_pres,local_profile, local_month)



outfile = "clim." + modelvarname + ".nc"
import scipy.io.netcdf as NC
ncOUT = NC.netcdf_file(outfile,"w")
ncOUT.createDimension("seas", nSeas)
ncOUT.createDimension("nSub", nSub)
ncOUT.createDimension("layers", nLayers)
ncOUT.createDimension("Statistics", nStat)

ncvar = ncOUT.createVariable("CLIMATOLOGY", "f", ('seas','nSub','layers','Statistics'))
ncvar[:] = CLIM
ncOUT.close()


import pylab as pl
PresCEN = [(l.bottom+l.top)/2  for l in LayerList]
outdir = "FIG_CLIM2/"

for it, seasreq in enumerate(REQUESTORS_LIST):
    for isub, sub in enumerate(basV2.P):
        
        Data_for_plot = CLIM[it,isub,:,:] #nStat x nLayers
        if np.all(np.isnan(Data_for_plot)) :
            print "No data for ", sub.name, seasreq.string
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

