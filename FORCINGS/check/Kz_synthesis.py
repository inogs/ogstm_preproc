import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''

    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'metrics_2d dir')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc')

    return parser.parse_args()

args = argument()


import numpy as np
from commons.mask import Mask
from commons.submask import SubMask
from commons.dataextractor import DataExtractor
from commons.Timelist import TimeList
from basins import V2 as OGS
from commons import season
from commons import timerequestors
from commons.utils import addsep




SeasonObj = season.season()
INPUTDIR=addsep(args.inputdir)
OUTDIR = addsep(args.outdir)


TheMask=Mask(args.maskfile)
jpk, jpj, jpi = TheMask.shape
TheMask.cut_at_level(0)

TL = TimeList.fromfilenames(None, INPUTDIR, "metrics*nc",prefix="metrics.")
dtype = [(sub.name, bool) for sub in OGS.P]
SUBlist=[ sub.name for sub in OGS.P.basin_list ]
nSub = len(SUBlist)


SUB = np.zeros((jpj,jpi),dtype=dtype)
SUBPoints = np.zeros((nSub,), dtype=int)

for isub, sub in enumerate(OGS.P):
    if sub.name=="med":continue
    m = SubMask(sub,maskobject=TheMask)
    SUB[sub.name] = m.mask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]
    SUBPoints[isub] = SUB[sub.name].sum()

ind=OGS.P.basin_list.index(OGS.med)
SUBPoints[ind] = SUB[OGS.med.name].sum()


def store_sub_arrays(TL, indexes, varname):
    
    '''
    returns
    a List of ndarrays, one for every subbasin
    Every array contains the values of varname, of its subbasin, queued in time
    '''
    nFrames = len(indexes)
    LIST_of_ndarrays=[]
    for isub, sub in enumerate(SUBlist):
        LIST_of_ndarrays.append(np.zeros((SUBPoints[isub]*nFrames,),np.float32))
    
    #Timeseries = np.zeros((nFrames_season,jpj,jpi),np.float32)
    for iframe, ind in enumerate(indexes):
        filename=TL.filelist[ind]
        M2d=DataExtractor(TheMask,filename,varname).values
        #Timeseries[iframe,:] = r
        for isub, sub in enumerate(OGS.P):
            l = LIST_of_ndarrays[isub]
            nPoints=SUBPoints[isub]
            mask = SUB[sub.name]
            l[nPoints*iframe:nPoints*(iframe+1)]  = M2d[mask]
    return LIST_of_ndarrays
def get_percentiles(LIST_of_arrays): 
    '''
    returns a 2d array [nPerc, nSub] of percentiles
    '''
    perc=[1,5,10,25,50,75,90,95,99]
    nperc=len(perc)
    PERC = np.zeros((nperc,nSub),np.float32)
    for isub, sub in enumerate(SUBlist):
        if (SUBPoints[isub] > 0 ):
            l = LIST_of_arrays[isub]
            good = ~np.isnan(l)
            PERC[:,isub] = np.percentile(l[good],perc)
    return PERC 

for iseas in range(4):
        
    req = timerequestors.Season_req(2019,iseas,SeasonObj)
    print(req)
    indexes, w = TL.select(req)
    LIST_of_ndarrays  = store_sub_arrays(TL, indexes, "Ved_150")
    PERC_surf = get_percentiles(LIST_of_ndarrays)
    outfile="%sPercentiles.%s.surf.npy" %(OUTDIR,iseas)
    np.save(outfile,PERC_surf)
    
    LIST_of_ndarrays  = store_sub_arrays(TL, indexes, "Ved_500")
    PERC_deep = get_percentiles(LIST_of_ndarrays)
    outfile="%sPercentiles.%s.deep.npy" %(OUTDIR,iseas)
    np.save(outfile,PERC_deep)

        
            
            
        
         
    