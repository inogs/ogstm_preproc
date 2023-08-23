import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates histograms

    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/gpfs/work/OGS_prod_0/OPA/V5C/devel/wrkdir/2/MODEL/FORCINGS/')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc')
    parser.add_argument(   '--top', '-t',
                                type = str,
                                required = True,
                                help = 'Top of the layer, in meters')
    parser.add_argument(   '--bottom', '-b',
                                type = str,
                                required = True,
                                help = 'Bottom of the layer, in meters')
    

    return parser.parse_args()

args = argument()

from commons.mask import Mask
from commons.dataextractor import DataExtractor
import numpy as np
from commons.Timelist import TimeList
from commons.utils import addsep
from commons.layer import Layer
from commons import netcdf4
import sys
from commons.submask import SubMask
from basins import V2 as OGS

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outdir)
TheMask=Mask(args.maskfile)


jpk, jpj, jpi = TheMask.shape

TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", prefix="U")



layer=Layer(args.top,args.bottom)
iz_top=TheMask.getDepthIndex(layer.top)
iz_bot=TheMask.getDepthIndex(layer.bottom)

SUBlist=[ sub.name for sub in OGS.P.basin_list ]
nSub = len(SUBlist)
mydtype = [(sub.name, bool) for sub in OGS.P.basin_list ]
SUBM = np.zeros((iz_bot-iz_top,jpj,jpi),dtype=mydtype)
SUBPoints = np.zeros((nSub,), dtype=int)

index_med=SUBlist.index('med')
for sub in SUBlist:
    index= SUBlist.index(sub)
    if index==index_med: continue
    basin = OGS.P.basin_list[index]
    s=SubMask(basin,maskobject = TheMask)
    SUBM[sub] = s.mask[iz_top:iz_bot,:,:]
    SUBPoints[index] = SUBM[sub].sum()
#  SUBM med is calculated as OR of all the others

for sub in OGS.Pred.basin_list[:-1]: # removing atlantic
    SUBM['med']=(SUBM['med'] | SUBM[sub.name])
SUBPoints[index_med] = SUBM['med'].sum()



bins=[0,1.e-7,1.e-6,1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1, 10]
nbins=len(bins)

from commons import season
from commons import timerequestors
SeasonObj = season.season()

for iseas in range(4):
    req = timerequestors.Season_req(2019,iseas,SeasonObj)
    outfile = "%sHistograms.%s.%s.npy" %(OUTPUTDIR, req.longname,layer.string())
    print(outfile)
    
     
    indexes, w = TL.select(req)
    nFrames = len(indexes)
    LIST_of_ndarrays=[]
    for isub, sub in enumerate(SUBlist):
        LIST_of_ndarrays.append(np.zeros((SUBPoints[isub]*nFrames,),np.float32))  
    
    for iframe,ind in enumerate(indexes):
        timestr = TL.Timelist[ind].strftime("%Y%m%d-%H:%M:%S")
        print(timestr,flush=True)
        filenameW=INPUTDIR + "W" + timestr + ".nc"
        K=DataExtractor(TheMask,filenameW,"votkeavt").values[iz_top:iz_bot,:,:]
        for isub, sub in enumerate(SUBlist):
            l = LIST_of_ndarrays[isub]
            nPoints = SUBPoints[isub]
            mask = SUBM[sub]
            l[nPoints*iframe:nPoints*(iframe+1)] = K[mask]
    
    HIST = np.zeros((nbins-1,nSub),int)
    for isub, sub in enumerate(SUBlist):
        HIST[:,isub],_ = np.histogram(LIST_of_ndarrays[isub],bins)
    np.save(outfile,HIST)
    
    
    


