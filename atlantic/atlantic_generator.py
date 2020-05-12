from commons.mask import Mask
import numpy as np
from commons import netcdf4

INPUTDIR="/gpfs/scratch/userexternal/vdibiagi/ogstm_boundary_conditions/atlantic/gen_inputs_BC/PROFILES/"
OUTPUTDIR="/gpfs/scratch/userexternal/vdibiagi/ogstm_boundary_conditions/atlantic/BC_atl_boundary/"

#check 
OpenMask=Mask('/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc')

jpk, jpj, jpi = OpenMask.shape
#check 
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
    outfile=OUTPUTDIR + "OPE_yyyy0630-00:00:00.nc"
    netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)


