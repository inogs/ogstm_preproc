# similar to aveScan
import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.basins import V2 as OGS
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.Timelist import TimeList
import sys
import pylab as pl

TheMask=Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask_CMCC.nc')
INPUTDIR="/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/output/"
OUTDIR = "/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/Ved_ratio/"

jpk, jpj,jpi = TheMask.shape
bool_dtype=[(sub.name,bool)    for sub in OGS.P]
Mask_0 = TheMask.cut_at_level(0)
SUBM = np.zeros((jpj, jpi), dtype=bool_dtype)
for isub, sub in enumerate(OGS.Pred):
    S=SubMask(sub, maskobject=Mask_0)
    SUBM[sub.name]=S.mask
    SUBM['med']=(SUBM['med'] | SUBM[sub.name])
    
TL = TimeList.fromfilenames(None, INPUTDIR, "metrics*nc",prefix="metrics.")
nFrames = TL.nTimes
nSub=len(OGS.P.basin_list)

varname="Ved_ratio"

Mratio = np.zeros((nFrames,nSub),dtype=np.float32)
Mcounter = np.zeros((nFrames,nSub),dtype=np.float32)


for iframe, filename in enumerate(TL.filelist):
    print(filename)
    r=DataExtractor(TheMask,filename,varname).values
    THRESHOLD = 5
    for isub, sub in enumerate(OGS.P):
        mask=SUBM[sub.name]
        ii = mask & (~np.isnan(r))
        nPoints = ii.sum()
        bad = r[ii]> THRESHOLD
        nBad = bad.sum()
        #M[iframe,isub] = np.nanmean(r[mask]) # because is base on deep layers
        Mratio[iframe,isub] = nBad/nPoints

    THRESHOLD = 2
    for isub, sub in enumerate(OGS.P):
        mask=SUBM[sub.name]
        nPoints = mask.sum()
        bad = r[mask]> THRESHOLD
        nBad = bad.sum()
        Mcounter[iframe,isub] = nBad/nPoints




for isub, sub in enumerate(OGS.P):
    fig,ax = pl.subplots()
    ax.plot(TL.Timelist,Mratio[:,isub])
    fig.suptitle(sub.name)
    outfile = OUTDIR + sub.name + ".png"
    fig.savefig(outfile)
    pl.close(fig)
    
OUTDIR = "/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/Anom_counter/"

       

for isub, sub in enumerate(OGS.P):
    fig,ax = pl.subplots()
    ax.plot(TL.Timelist,Mcounter[:,isub])
    fig.suptitle(sub.name)
    outfile = OUTDIR + sub.name + ".png"
    fig.savefig(outfile)
    pl.close(fig)


from bitsea.commons.utils import writetable
from bitsea.commons import timerequestors
from bitsea.commons import season
SUBLIST=[ sub.name for sub in OGS.P ]
SeasonObj=season.season()

nSUB = len(SUBLIST)
CounterTable=np.zeros((nSub,3),np.float32)
CounterTable[:,0] = Mcounter[:,:].mean(axis=0)

#writetable(OUTDIR + "annual.txt",Mcounter[:,:].mean(axis=0).reshape(18,1),SUBLIST,['mean'])

req=timerequestors.Season_req(2019,0,SeasonObj)
ii,w = TL.select(req)
CounterTable[:,1] = Mcounter[ii,:].mean(axis=0)

#writetable(OUTDIR + "winter.txt",Mcounter[ii,:].mean(axis=0).reshape(18,1),SUBLIST,['mean'])

req=timerequestors.Season_req(2019,2,SeasonObj)
ii,w = TL.select(req)
CounterTable[:,2] = Mcounter[ii,:].mean(axis=0)

writetable(OUTDIR + "mean.txt",CounterTable,SUBLIST,['annual','winter','summer'])



    
