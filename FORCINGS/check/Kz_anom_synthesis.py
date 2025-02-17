import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Kz_Anomalies.txt and Kz_Anomalies.png, about Vertical Eddy Diffusivity
    The two file contains is a table [sub, season ] of ratio
    column with anomalous values / total of watercolumns
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
                                help = 'meshmask_INGV.nc')

    return parser.parse_args()

args = argument()



import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.Timelist import TimeList
from bitsea.basins import V2 as OGS
from bitsea.commons import season
from bitsea.commons import timerequestors
from bitsea.commons.utils import writetable
from bitsea.commons.utils import addsep
import matplotlib
matplotlib.use('Agg')
import pylab as pl

SeasonObj = season.season()


TheMask=Mask.from_file(args.maskfile)
jpk, jpj, jpi = TheMask.shape
TheMask.cut_at_level(0)

dtype = [(sub.name, bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)

for isub, sub in enumerate(OGS.Pred):
    m = SubMask(sub, TheMask)
    SUB[sub.name] = m.mask[0,:,:]
    SUB['med'] = SUB['med'] | SUB[sub.name]


INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outdir)


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

writetable(OUTPUTDIR + 'Kz_Anomalies.txt', Mout, rows_names_list, column_names, fmt="%5.3f\t")


pl.close('all')
fig,ax=pl.subplots()
fig.set_size_inches(5,5)
ax.set_position([0.08, 0.13, 0.78, 0.78])
cmap=pl.get_cmap('viridis',10)
im=ax.imshow(Mout*100, aspect='auto',extent=[-.5,3.5,-0.5,17.5], cmap=cmap)

im.set_clim(0, 20)

cbar_ticks_list = np.arange(0,22,2)
cbar = fig.colorbar(im,ticks=cbar_ticks_list)

ax.set_yticks(np.arange(nSub))
ax.set_yticklabels(rows_names_list[-1::-1])
ax.set_xticks(np.arange(4))
ax.set_xticklabels(column_names)

for isub,sub in enumerate(OGS.P):
    for iseas in range(4):
        y=nSub-1-isub
        x=iseas
        value=Mout[isub,iseas]*100
        color="w"
        if value>12: color="k"
        text="%3.2f" %(value)
        ax.text(x,y,text,color=color, fontsize=8, ha='center', va='center')
fig.suptitle("% Kz anomalies in 150 - 500m")
fig.savefig(OUTPUTDIR + 'Kz_anomalies.png',dpi=150)



