from commons.mask import Mask
import numpy as np
from commons import netcdf4
from dardanelles.analysis_config import VARLIST, mydtype

#BOUNDARY_CONCENTRATION =  np.load('../dard_boundary_values.npy')
BOUNDARY_CONCENTRATION= np.zeros((1,),dtype=mydtype)
BOUNDARY_CONCENTRATION['N5s'] = 2.0
BOUNDARY_CONCENTRATION['N1p'] = 0.065 # mmol/m3
BOUNDARY_CONCENTRATION['N3n']=   1.3 # mol/m3
BOUNDARY_CONCENTRATION['O3c']= 28700 # mg/m3
BOUNDARY_CONCENTRATION['O3h']=  2800 # mmol/m3

def dump_files(OpenMask, OUTPUTDIR):

    jpk, jpj, jpi = OpenMask.shape
    I=841

    for var in VARLIST:

        BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
        for k in range(jpk):
            BOUNDARY[k,235:238,I] = BOUNDARY_CONCENTRATION[var][0]
        BOUNDARY[~OpenMask.mask]=1.e+20
        outfile=OUTPUTDIR + "OPE_yyyy0630-00:00:00.nc"
        netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)

    
             
        
    