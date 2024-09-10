import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''Read datasets and compute statistics 
                               over vertical layers in the 16 Med subbasins 
                              
                               ''')
    parser.add_argument(   '--parentoutput', '-p',
                                type = str,
                                required = True,
                                help ='/gpfs/scratch/userexternal/vdibiagi/gen_inputs_IC/fromDataset/'
                                )
    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                help ='N1p'
                                )
  
    parser.add_argument(   '--start_time', '-s',
                                type = str,
                                required = True,
                                help ='1997'
                                )
    parser.add_argument(   '--end_time', '-e',
                                type = str,
                                required = True,
                                help ='2015'
                                )

        
    return parser.parse_args()


args = argument()

from bitsea.commons.time_interval import TimeInterval
from bitsea.basins import V2 as OGS
import doStats as doStats
from bitsea.commons.layer import Layer
from bitsea.commons import timerequestors
from bitsea.commons import season
import numpy as np
import matplotlib.pyplot as plt
from bitsea.instruments.var_conversions import NUTRVARS
from bitsea.commons.utils import addsep
import os

varDataName = args.varname
firstYear = args.start_time
lastYear =args.end_time
parentOutput  = addsep(args.parentoutput)
os.system("mkdir -p " + parentOutput)
OUTPUTDIR = addsep(parentOutput + "/" + firstYear + "_" + lastYear) 
os.system("mkdir -p " + OUTPUTDIR)

firstYearExcluded = str(int(lastYear)+1)
TI=TimeInterval(firstYear + '0101',firstYearExcluded + '0101','%Y%m%d')

#Season selection
S=season.season()

varName, DataSet = doStats.DatasetInfo(varDataName)

LayerList=[Layer(0,25),Layer(25,50), Layer(50,75), Layer(75,100), Layer(100,125), 
Layer(125,150), Layer(150, 200), Layer(200,400), Layer(400,600), Layer(600, 800),
Layer(800, 1000), Layer(1000, 1500), Layer(1500, 2000), Layer(2000, 2500)]

print "Computation on subbasins..."

MEAN, STD, NVALS, NPROF_SEAS= doStats.get_yearsStats(varDataName, OGS.P.basin_list, LayerList, TI)  
# MEAN, STD have dimensions: nSub, nLayer, 5 (one annual profile and 4 seasonal profiles)

# midpoint of the layers 
nLayers = len(LayerList)
depthLayer = np.zeros(nLayers, np.float32)*np.nan
for ilayer, layer in enumerate(LayerList):
    depthLayer[ilayer] = 0.5*(layer.top+layer.bottom)


# 0:annual, 1:winter, 2:spring, 3:summer, 4:fall
seasNameOut=[]
for iseas in range(0,5): 
	MMEAN =np.zeros((14,),np.float32)
	MSTD =np.zeros((14,),np.float32)
	MNVALS=np.zeros((14,),np.float32)

	if iseas==0:
		seasNameOut.append('annual')
	else:
		reqSeas=timerequestors.Clim_season(iseas-1,S) 
		seasNameOut.append(reqSeas.seasonobj.SEASON_LIST_NAME[iseas-1])
	print "season name: ", seasNameOut
	for isub in range(0,17): 
		#print OGS.P.basin_list[isub].name
		mean=MEAN[isub,:,iseas] 
		std=STD[isub,:,iseas]
		nvals=NVALS[isub,:,iseas]

		if isub>0:
			MMEAN = np.vstack((MMEAN,mean))
			MSTD = np.vstack((MSTD,std))
			MNVALS = np.vstack((MNVALS,nvals))
		else:
			MMEAN=mean
			MSTD=std
			MNVALS=nvals
	np.savetxt(OUTPUTDIR+seasNameOut[iseas]+'_vertprof_'+varDataName+'.txt',np.concatenate((MMEAN,MSTD),axis=1))
	np.savetxt(OUTPUTDIR+seasNameOut[iseas]+'_nvals_'+varDataName+'.txt',MNVALS)
	np.savetxt(OUTPUTDIR+'vertprof_depthLayer.txt',depthLayer)
	np.savetxt(OUTPUTDIR+'numberOfProfiles_'+varDataName+'.txt',NPROF_SEAS,fmt='%5.0f')

