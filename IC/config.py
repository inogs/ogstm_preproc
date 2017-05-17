from basins.basin import ComposedBasin
from basins import V2
from commons.layer import Layer
import numpy as np
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