from bitsea.commons.mask import Mask
import numpy as np
from bitsea.commons import netcdf4

maskfile="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/meshmask.nc"
TheMask=Mask(maskfile)


I=841
J=[235,236,237]

jpk,jpj,jpi=TheMask.shape

M=np.zeros((jpk,jpj,jpi),np.float32)

M[:,235:238,I] = 1


netcdf4.write_3d_file(M, 'geom', "OPEN_geom.nc", TheMask, compression=True)