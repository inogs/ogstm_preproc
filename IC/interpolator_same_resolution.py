import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Interpolates with nearest algorithm restarts between very similar meshes
    Meshes are supposed to have the same definition and be different from some land/sea points
 
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'original restarts')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output directory of generated restarts'
                                )
    parser.add_argument(   '--origmask',
                                type = str,
                                required = False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V4/meshmask.nc",
                                help = 'Path of the mask file'
                                )
    parser.add_argument(   '--newmask',
                                type = str,
                                required = False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V4.1/meshmask.nc",
                                help = 'Path of the mask file'
                                )    

    return parser.parse_args()

args = argument()


import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
import os,glob
import scipy.interpolate
import netCDF4
from IC import RSTwriter
from bitsea.commons.utils import addsep

TheMask_new=Mask.from_file(args.newmask)
TheMask_orig =Mask.from_file(args.origmask)
ORIGDIR= addsep(args.inputdir) #"/gpfs/scratch/userexternal/ateruzzi/DA_FloatNut/RUN_REFnew/wrkdir/MODEL/RESTARTS/"
OUTDIR = addsep(args.outdir) #"/gpfs/scratch/userexternal/ateruzzi/DA_FloatNut/RUN_REFnew/INIT/"

jpk, jpj,jpi = TheMask_new.shape
filelist=glob.glob(ORIGDIR + "RST.20150101*nc")

for filename in filelist:
    basename=os.path.basename(filename)
    outfile = OUTDIR + basename
    _,_,var,_=basename.rsplit(".")
    RST_old = DataExtractor(TheMask_orig,filename, "TRN"+var).values
    RST_old[~TheMask_orig.mask]=np.nan
    RST=RST_old[:jpk,:,:]
    
    for k in range(jpk):
        rst=RST[k,:,:]
        goods  = TheMask_orig.mask[k,:,:]
        Jgoods, Igoods = np.nonzero(goods)
        nP = len(Jgoods)
        points = np.zeros((nP,2),dtype=np.float32)
        points[:,0] = Jgoods
        points[:,1] = Igoods
        values = rst[goods]
        
        tofill = (~TheMask_orig.mask[k,:,:] & TheMask_new.mask[k,:,:])
        J,I = np.nonzero(tofill)
        nP = len(J)
        xi = np.zeros((nP,2),dtype=np.float32)
        xi[:,0] = J
        xi[:,1] = I
        
        V= scipy.interpolate.griddata(points, values, xi, "nearest")
        rst[tofill] = V
    RST[~TheMask_new.mask] = 1.e+20
    check = np.isnan(RST[TheMask_new.mask])
    print var + " : number of nans: ", check.sum()
    RSTwriter(outfile, var, RST, TheMask_new)
        