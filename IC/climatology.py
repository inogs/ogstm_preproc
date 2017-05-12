from basins.basin import ComposedBasin
from basins import V2
from commons.layer import Layer
from commons.time_interval import TimeInterval
from static.Nutrients_reader import NutrientsReader
from static.Carbon_reader import CarbonReader
from instruments.var_conversions import NUTRVARS

import numpy as np

N=NutrientsReader()
C=CarbonReader()

def setvarname(var):
    if var in ['N1p','N3n','O2o','N5s']:
        return NUTRVARS[var]
    if var == 'O3h': return 'DIC'
    if var == 'O3c': return 'ALK'


DATASET_DICT = {'N1p':N, 'N3n':N, 'O2o':N, 'N5s':N, 'O3c':C, 'O3h':C}

VARLIST=['N1p','N3n','O2o','N5s','O3h','O3c']

TI = TimeInterval("1950","2050","%Y")


PresDOWN=np.array([25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[]
top = 0
for bottom in PresDOWN:
    LayerList.append(Layer(top, bottom))
    top = bottom

Atl = ComposedBasin('atl',[V2.alb,   V2.atl],'Alboran & Atlantic')
NWM = ComposedBasin('nwm',[V2.nwm,  V2.tyr1],'North West Med')
SWM = ComposedBasin('swm',[V2.tyr2, V2.swm2],'South West Med')
AEG = ComposedBasin('aeg',[V2.aeg,  V2.lev1],'Aegean Sea')
ION = ComposedBasin('ion',[V2.ion1, V2.ion2],'Ionian Sea')
ICdef = ComposedBasin('ICdef',[Atl, NWM, SWM, V2.swm1, V2.adr, AEG, ION, V2.ion3, V2.lev],'Gruped Subbasin for Climatology for Restarts')




nLayers = len(LayerList)
nVars   = len(VARLIST)
CLIM = np.zeros((nVars,nLayers), np.float32)*np.nan

for ivar, modelvarname in VARLIST:
    Dataset =  DATASET_DICT[modelvarname]
    var     =  setvarname(modelvarname)
    for isub, sub in enumerate(ICdef):
        print " --- * ", sub.name
        Profilelist =Dataset.Selector(var, TI, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        for p in Profilelist: 
            pres, profile, Qc = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<=layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[ivar, ilayer] = Values[ii].mean()



