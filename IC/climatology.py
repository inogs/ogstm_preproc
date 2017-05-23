from config import ICdef, LayerList
from commons.time_interval import TimeInterval
from static.Nutrients_reader import NutrientsReader
from static.Carbon_reader import CarbonReader
from instruments.var_conversions import NUTRVARS

import numpy as np

N=NutrientsReader()
C=CarbonReader()



def DatasetInfo(modelvarname):
    '''
    Argument:
    * modelvarname * string, like 'N1p'
    Returns:
    * var     * variable name to access dataset
    * dataset * NutrientsReader or Carbonreader object
    '''
    if modelvarname in ['N1p','N3n','O2o','N5s']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    if modelvarname == 'O3h':
        var = 'ALK'
        dataset = C
    if modelvarname == 'O3c':
        var ='DIC'
        dataset = C
    return var, dataset



TI = TimeInterval("1950","2050","%Y")

nLayers = len(LayerList)
nSub    = len(ICdef.basin_list)


def get_climatology(modelvarname):
    CLIM = np.zeros((nSub, nLayers), np.float32)*np.nan
    var, Dataset = DatasetInfo(modelvarname)
    for isub, sub in enumerate(ICdef):
        Profilelist =Dataset.Selector(var, TI, sub)
        Pres  =np.zeros((0,),np.float32)
        Values=np.zeros((0,),np.float32)
        for p in Profilelist: 
            pres, profile, _ = p.read(var)
            Pres   = np.concatenate((Pres,pres))
            Values = np.concatenate((Values,profile))
        for ilayer, layer in enumerate(LayerList):
            ii = (Pres>=layer.top) & (Pres<=layer.bottom)
            if (ii.sum()> 1 ) :
                CLIM[isub, ilayer] = Values[ii].mean()
    return CLIM
