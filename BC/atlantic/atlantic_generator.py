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

from bitsea.commons.mask import Mask
import numpy as np
from bitsea.commons import netcdf4
from bitsea.commons.utils import addsep
import os, glob
from smoother_1d import smoother_1d

INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR  = addsep(args.outputdir)
RESTARTDIR = addsep(args.rst)
MASKFILE=args.maskfile

os.system("mkdir -p " + OUTPUTDIR)

OpenMask=Mask.from_file(MASKFILE)
nav_lon=OpenMask.xlevels
B=OpenMask.bathymetry_in_cells()
nav_lev=OpenMask.zlevels
jpk, jpj, jpi = OpenMask.shape
I=1

lon_max_Atl=-6.125
diff=abs(nav_lon[0,:]-lon_max_Atl)
idx= np.where(diff == np.min(diff))
jpi_Gib=idx[0][0] 
arrAtl=B[:,0:jpi_Gib]
maxLevel = np.max(arrAtl) 
position = np.where(arrAtl == np.max(arrAtl))
jpj_p=position[0][0]
jpi_p=position[1][0]

os.chdir(RESTARTDIR)
PATH_NAME='RST.*.nc'
fileLIST = glob.glob(PATH_NAME)
ALLVARS=[f[-6:-3] for f in fileLIST]
mydtype=[(var,np.float32) for var in ALLVARS]
BOUNDARY_CONCENTRATION= np.zeros((jpk,),dtype=mydtype)

VARLIST=['N1p', 'N3n', 'N5s','O2o','O3c','O3h','Hg2','Hg0','MHg','DHg','P1h','P2h','P3h','P4h','Z6h','Z5h','Z4h','Z3h']

for var in ALLVARS:
    if var in VARLIST:
        filename=INPUTDIR + var + ".nc"
        profile=netcdf4.readfile(filename,var)
        BOUNDARY_CONCENTRATION[var][:]=profile
    else:
        filename=RESTARTDIR + fileLIST[0][0:-6] + var + ".nc"
        profile=netcdf4.readfile(filename,"TRN"+var)
        profile_sel=profile[0][:,jpj_p,jpi_p] # profile in the deepest Atlantic point
        profile_smoothed=profile_sel[0:maxLevel]
        for itime in range(2):
            profile_smoothed=smoother_1d(profile_smoothed,nav_lev[0:maxLevel])
        check=np.max(profile_smoothed)
        BOUNDARY_CONCENTRATION[var][0:maxLevel]=profile_smoothed
        BOUNDARY_CONCENTRATION[var][maxLevel:-1]=profile_smoothed[-1]
    check=np.max(BOUNDARY_CONCENTRATION[var][:])        
    print(var, "max value: ", check)
    BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
    for k in range(jpk):
        BOUNDARY[k,:,I] = BOUNDARY_CONCENTRATION[var][k]
    BOUNDARY[~OpenMask.mask]=1.e+20
    outfile=OUTPUTDIR + "ATL_yyyy0630-00:00:00.nc"
    netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)


