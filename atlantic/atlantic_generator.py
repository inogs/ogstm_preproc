import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''Creates 3D files of BC in Atlantic,
                               starting from netCDF profiles
                              
                               ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/inputs'
                                )
    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/OUTPUT'
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc'
                                )

    return parser.parse_args()


args = argument()

from commons.mask import Mask
import numpy as np
from commons import netcdf4
from commons.utils import addsep
import os

INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR  = addsep(args.outputdir)
MASKFILE=args.maskfile

os.system("mkdir -p " + OUTPUTDIR)

OpenMask=Mask(MASKFILE)

jpk, jpj, jpi = OpenMask.shape
I=1

VARLIST=['N1p', 'N3n', 'N5s','O2o','O3c','O3h']
mydtype=[(var,np.float32) for var in VARLIST]

BOUNDARY_CONCENTRATION= np.zeros((jpk,),dtype=mydtype)

for var in VARLIST:

    filename=INPUTDIR + var + ".nc"
    profile=netcdf4.readfile(filename,var) 
    BOUNDARY_CONCENTRATION[var][:]=profile
    BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
    for k in range(jpk):
        BOUNDARY[k,:,I] = BOUNDARY_CONCENTRATION[var][k]
    BOUNDARY[~OpenMask.mask]=1.e+20
    outfile=OUTPUTDIR + "ATL_yyyy0630-00:00:00.nc"
    netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)


