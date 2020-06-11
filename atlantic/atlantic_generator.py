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
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/OUTPUT/ALLVARS'
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc'
                                )
    parser.add_argument(   '--rst', '-r',
                                type = str,
                                required = True,
                                help ='/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/IC/RST'
                                )

    return parser.parse_args()


args = argument()

from commons.mask import Mask
import numpy as np
from commons import netcdf4
from commons.utils import addsep
import os, glob
from scipy import interpolate
from smoother_1d import smoother_1d

INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR  = addsep(args.outputdir)
RESTARTDIR = addsep(args.rst)
MASKFILE=args.maskfile

os.system("mkdir -p " + OUTPUTDIR)

OpenMask=Mask(MASKFILE)
nav_lev=OpenMask.zlevels
jpk, jpj, jpi = OpenMask.shape
I=1

os.chdir(RESTARTDIR)
PATH_NAME='RST.*.nc'
fileLIST = glob.glob(PATH_NAME)

VARLIST=['N1p', 'N3n', 'N5s','O2o','O3c','O3h']

ALLVARS=[f[-6:-3] for f in fileLIST]
mydtype=[(var,np.float32) for var in ALLVARS]
BOUNDARY_CONCENTRATION= np.zeros((jpk,),dtype=mydtype)
for var in ALLVARS:
	if var in VARLIST:
    		filename=INPUTDIR + var + ".nc"
    		profile=netcdf4.readfile(filename,var) 
    		BOUNDARY_CONCENTRATION[var][:]=profile
        else:
                filename=RESTARTDIR + fileLIST[0][0:-6] + var + ".nc"
                profile=netcdf4.readfile(filename,"TRN"+var)
                profile_sel=profile[0][:,102,2] #atlantic point
                profile_smoothed=profile_sel
                for itime in range(2):
			profile_smoothed=smoother_1d(profile_smoothed,nav_lev)
                BOUNDARY_CONCENTRATION[var][:]=profile_smoothed
    	BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
    	for k in range(jpk):
        	BOUNDARY[k,:,I] = BOUNDARY_CONCENTRATION[var][k]
        BOUNDARY[~OpenMask.mask]=1.e+20
	outfile=OUTPUTDIR + "ATL_yyyy0630-00:00:00.nc"
	netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)
		

