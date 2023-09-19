import numpy as np
from commons.mask import Mask
from commons.submask import SubMask
from commons.dataextractor import DataExtractor
from commons.Timelist import TimeList
from basins import V2 as OGS
from commons import season
from commons import timerequestors
from commons.utils import writetable

SeasonObj = season.season()


TheMask=Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask_CMCC.nc')
jpk, jpj, jpi = TheMask.shape
TheMask.cut_at_level(0)

dtype = [(sub.name, bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)

for isub, sub in enumerate(OGS.Pred):
    m = SubMask(sub,maskobject=TheMask)
    SUB[sub.name] = m.mask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]


INPUTDIR="/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/output/"
OUTPUTDIR="/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/synthesis/"


TL = TimeList.fromfilenames(None, INPUTDIR, "metrics*nc",prefix="metrics.")
nSub=len(OGS.P.basin_list)

Mout = np.zeros((nSub,4),np.float32)
THRESHOLD = 2 

for iseas in range(4):
    M2d_sum = np.zeros((jpj,jpi),int)
    req = timerequestors.Season_req(2019,iseas,SeasonObj)
    print(req)
    indexes, w = TL.select(req)
    nFrames_season = len(indexes)

    for ind in indexes:
        filename=TL.filelist[ind]
        r=DataExtractor(TheMask,filename,'anom_counter').values
        Affected_columns = (r > THRESHOLD).astype(int)
        M2d_sum += Affected_columns
    for isub, sub in enumerate(OGS.P):
        ii = SUB[sub.name]
        nColumnsSub = ii.sum()
        Mout[isub,iseas] = M2d_sum[ii].sum()/(nFrames_season * nColumnsSub) 

rows_names_list=[sub.name for sub in OGS.P]
column_names = ['win','spr','sum','fal']

writetable(OUTPUTDIR + 'Anomalies.txt', Mout, rows_names_list, column_names, fmt="%5.3f\t")
    
    

    
