import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''

    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    return parser.parse_args()

args = argument()

import numpy as np
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
from bitsea.basins import V2 as basV2
from bitsea.static import climatology
from bitsea.commons.utils import addsep
from bitsea.commons import timerequestors
import xarray as xr

LayerList_metrics = [Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]


PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList_plot=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

OUTDIR = addsep(args.outdir)
TI = TimeInterval("19990101","20240101","%Y%m%d")

TheMask = Mask.from_file(args.maskfile)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels


Req=timerequestors.Generic_req(TI)

VARLIST1=['O2o','ALK','DIC','pH','N4n','pCO2']
VARLIST2=['N1p','N3n','N5s']
SUBlist = basV2.Pred.basin_list[:]
# remove Atlantic buffer from the list
SUBlist.remove(SUBlist[-1])
nSub = len(SUBlist)
print(nSub)

# Crea il nuovo dataset
def dumpfile(filename,z_levels, Mean,Std,included_subds, nprofiles=None):

    new_ds = xr.Dataset(
            {
                "mean": (("nsubds", "depth"), Mean),
                "std": (("nsubds", "depth"), Std),
                "total_profs": (("nsubds",), nprofiles),
                "included_profs_flag":(("nsubds",), included_subds),
            },
            coords={
                "nsubds": [sub.name for sub in SUBlist],
                "depth": z_levels,
            },
            attrs={
                "description": "Annual climatology in open sea from MedBGCins in 1999-2024"
            },
        
        )
    new_ds.to_netcdf(filename,format="NETCDF4")


#Plot section (QC=false)
z_clim_plot = np.array([-(el.bottom+el.top)/2  for el in LayerList_plot])
nLayers = len(LayerList_plot)

included=np.ones((nSub,),int)

for ivar, var in enumerate(VARLIST1):
    CLIM_REF_static,STD_Ref_static,_,Nprofs = climatology.get_climatology_open(var,SUBlist, LayerList_plot, TheMask, useLogistic=False, basin_expand=False,QC=False, climatology_interval=TI)
    outfile =OUTDIR + var + "_clim_plot.nc"
    dumpfile(outfile,z_clim_plot, CLIM_REF_static, STD_Ref_static,included, Nprofs)

for ivar, var in enumerate(VARLIST2):
    if var=='N1p':
        CLIM_REF_static,STD_Ref_static,_,Nprofs = climatology.get_climatology_open(var, SUBlist, LayerList_plot, TheMask, useLogistic=True,startpt=np.asarray([0.1, 0.1, 1000, 0.4],dtype=np.float64), basin_expand=False, QC=False, climatology_interval=TI)
    elif var=='N3n' or var=='N5s':
        CLIM_REF_static,STD_Ref_static,_,Nprofs= climatology.get_climatology_open(var, SUBlist, LayerList_plot, TheMask, useLogistic=True,startpt=np.asarray([0.1, 0.1, 500, 4],dtype=np.float64),basin_expand=False, QC=False, climatology_interval=TI)
    else:
        print("VARLIST2 can include only N1p, N3n, N5s")
    outfile =OUTDIR + var + "_clim_plot.nc"
    dumpfile(outfile,z_clim_plot, CLIM_REF_static, STD_Ref_static,included, Nprofs)

# Metrics section (QC=True)

z_clim_metrics = np.array([-(el.bottom+el.top)/2  for el in LayerList_metrics])
nLayers = len(LayerList_metrics)
for ivar, var in enumerate(VARLIST1):
    included=np.ones((nSub,),int)
    for isub, sub in enumerate(SUBlist):
        included[isub] = climatology.QualityCheck(var, sub)
    print("included=",included)

    CLIM_REF_static,STD_Ref_static,_,Nprofs = climatology.get_climatology_open(var,SUBlist, LayerList_metrics, TheMask, useLogistic=False, basin_expand=False,QC=True, climatology_interval=TI)
    outfile =OUTDIR + var + "_clim_metrics.nc"
    dumpfile(outfile,z_clim_metrics, CLIM_REF_static, STD_Ref_static,included,Nprofs)

for ivar, var in enumerate(VARLIST2):
    included=np.ones((nSub,),int)
    for isub, sub in enumerate(SUBlist):
        included[isub] = climatology.QualityCheck(var, sub)
    print("included=",included)

for ivar, var in enumerate(VARLIST2):
    if var=='N1p':
        CLIM_REF_static,STD_Ref_static,_,Nprofs = climatology.get_climatology_open(var, SUBlist, LayerList_metrics, TheMask, useLogistic=True,startpt=np.asarray([0.1, 0.1, 1000, 0.4],dtype=np.float64), basin_expand=False, QC=True, climatology_interval=TI)
    elif var=='N3n' or var=='N5s':
        CLIM_REF_static,STD_Ref_static,_,Nprofs= climatology.get_climatology_open(var, SUBlist, LayerList_metrics, TheMask, useLogistic=True,startpt=np.asarray([0.1, 0.1, 500, 4],dtype=np.float64),basin_expand=False, QC=True, climatology_interval=TI)
    else:
        print("VARLIST2 can include only N1p, N3n, N5s")

    outfile =OUTDIR + var + "_clim_metrics.nc"
    dumpfile(outfile,z_clim_metrics, CLIM_REF_static, STD_Ref_static,included,Nprofs)



