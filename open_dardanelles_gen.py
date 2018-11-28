from commons.mask import Mask
from commons.submask import SubMask
import numpy as np
from commons import netcdf4
from commons.dataextractor import DataExtractor
from basins.region import Polygon
from basins.basin import SimplePolygonalBasin

p = Polygon([27.0483398,26.674804,24.82910,23.400878,22.456054,22.3022460,22.456054,23.07128,23.994140,27.35595], \
            [39.960280,40.797177,41.178653,41.029643,40.680638,39.943436,39.027718,38.272688,37.544577,37.61423])
local_sub  = SimplePolygonalBasin('aeg', p,'North Aegean')



OpenMask=Mask('/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/meshmask.nc')
maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc"
TheMask=Mask(maskfile)
mask_200= TheMask.mask_at_level(200.)
S = SubMask(local_sub,maskobject=TheMask)
jpk, jpj, jpi = TheMask.shape
I=841


INPUTDIR="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/INPUT_for_Dardanelles_profiles/"
OUTPUTDIR="/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/INPUTS_for_MODEL/"

VARLIST=['N1p', 'N3n', 'O2o', 'N5s','O3c','O3h']
for var in VARLIST:
    inputfile=INPUTDIR + "ave.20171002-12:00:00." + var + ".nc"
    V=DataExtractor(TheMask,inputfile,var).values
    Profile=np.zeros((jpk,),np.float32)*np.nan
    for k in range(jpk):
        mask=S.mask[k,:,:] & mask_200
        values =  V[k,:,:]
        if np.any(mask):
            Profile[k] = values[mask].mean()

    BOUNDARY = np.zeros((jpk,jpj,jpi), np.float32)
    for k in range(jpk):
        BOUNDARY[k,235:238,I] = Profile[k]
    BOUNDARY[~OpenMask.mask]=1.e+20
    outfile=OUTPUTDIR + "OPE.yyyy0630-00:00:00.nc"
    netcdf4.write_3d_file(BOUNDARY, var, outfile, OpenMask,compression=True)

    
             
        
    