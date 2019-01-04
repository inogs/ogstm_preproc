from commons.mask import Mask
import numpy as np
from commons import netcdf4
from analysis_config import VARLIST

BOUNDARY_CONCENTRATION =  np.load('../dard_boundary_values.npy')
BOUNDARY_CONCENTRATION['N1p'] = 0.05 # mmol/m3
BOUNDARY_CONCENTRATION['N3n']=   1.3 # mol/m3
BOUNDARY_CONCENTRATION['O3c']= 29200 # mg/m3
BOUNDARY_CONCENTRATION['O3h']=  2850 # mmol/m3
OUTPUTDIR="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/INPUTS_for_MODEL/2/"
OUTPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/BC/"

OpenMask=Mask('/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/meshmask.nc')

jpk, jpj, jpi = OpenMask.shape
I=841



for var in VARLIST:

    BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
    for k in range(jpk):
        BOUNDARY[k,235:238,I] = BOUNDARY_CONCENTRATION[var][0]
    BOUNDARY[~OpenMask.mask]=1.e+20
    outfile=OUTPUTDIR + "OPE.yyyy0630-00:00:00.nc"
    netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)

    
             
        
    