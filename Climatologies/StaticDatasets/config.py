# Definitions of climatology parameters: 
# Layers, subregions, seasons

import basins.V2 as basV2
from commons.layer import Layer
from commons import season, timerequestors
import numpy as np

PresDOWN=np.array([25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[]
top = 0
for bottom in PresDOWN:
    LayerList.append(Layer(top, bottom))
    top = bottom

SeasonObj = season.season()
AnnualObj = season.season()
AnnualObj.setseasons(["0101"], ["annual"])

nSeas   = SeasonObj.get_seasons_number()
REQUESTORS_LIST=[]
for iSeas in range(nSeas):
    REQUESTORS_LIST.append(timerequestors.Clim_season(iSeas, SeasonObj))
REQUESTORS_LIST.append(timerequestors.Clim_season(0,AnnualObj)) # adding annual