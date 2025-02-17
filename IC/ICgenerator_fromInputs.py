import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''Creates 3D files of IC,
                               starting from netCDF profiles in the 16 Med subbasins
                              
                               ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help ='/gpfs/scratch/userexternal/vdibiagi/toGiorgio/RAN24/ICs/MODEL_units/'
                                )
    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/IC/preproc/IC/RST/'
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc'
                                )

    parser.add_argument(   '--date', '-d',
                                type = str,
                                required = True,
                                help ='19990101'
                                )     
        
    return parser.parse_args()


args = argument()
from bitsea.commons.mask import Mask
import numpy as np
import pylab as pl
from bitsea.commons.submask import SubMask
from IC import RSTwriter, smoother
from bitsea.commons.utils import getcolor
from bitsea.commons import netcdf4
from bitsea.basins import V2 as OGS
from bitsea.commons.utils import addsep
import os

INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR  = addsep(args.outputdir)
DATE=args.date
MASKFILE=args.maskfile

os.system("mkdir -p " + OUTPUTDIR)

def get_climatology(varname):
    filename=INPUTDIR + varname + ".nc"
    CLIM_17=np.zeros((17,125), np.float32)
    CLIM_17[:16,:]=netcdf4.readfile(filename,varname) 
    CLIM_17[-1,:] = CLIM_17[0,:] # copy alb in atl
    return CLIM_17


TheMask = Mask.from_file(MASKFILE)

jpk, jpj, jpi = TheMask.shape

VARLIST=['N1p','N3n','O2o','N5s','O3h','O3c']
sublist17=[ sub for sub in OGS.Pred]
SUB17 = OGS.ComposedBasin('med',sublist17,'Mediterranean Sea with Atlantic buffer')

nSub = len(sublist17)
# first of all, plots
for varname in VARLIST:
    CLIM = get_climatology(varname)
    out_img = OUTPUTDIR + varname + ".png"
    fig, ax = pl.subplots()
    for isub, sub in enumerate(SUB17):
        p = CLIM[isub,:]
        ax.plot(p, TheMask.zlevels,color=getcolor(nSub, isub), label=sub.name)
    ax.invert_yaxis()
    ax.legend()
    fig.savefig(out_img)
    pl.close(fig)

# then 3D arrays generation ---------
for varname in VARLIST:
    CLIM = get_climatology(varname)
    outfile = OUTPUTDIR + "RST." + DATE + "-00:00:00." + varname + ".nc"
    print(outfile)
    RST = np.zeros((jpk,jpj,jpi),np.double)
    for isub, sub in enumerate(SUB17):
        p = CLIM[isub,:]
        S =SubMask(sub, TheMask)
        for k in range(jpk):
            submask = S.mask[k,:,:]
            V = RST[k,:,:]
            V[submask] =p[k]
    print("smoother")
    RST_s = smoother(TheMask, RST)
    print("writer")
    check = np.isnan(RST_s[TheMask.mask])
    print("number of nans: ", check.sum())
    RSTwriter(outfile, varname, RST_s, TheMask)

