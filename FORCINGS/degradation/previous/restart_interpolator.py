import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
import netCDF4
import glob,os

def weighted_mean(values,weights):
    return (values*weights).sum()/weights.sum()

def interp3d(M3dfine, Maskfine,MaskCoarse,I,J):
    '''
    Degradates a 3d field defined on cell center, using area weighted mean.

    Arguments:
    * M3dfine    * 3d numpy array
    * Maskfine   * Mask object
    * MaskCoarse * Mask object
    * I          * numpy array 1d of integers, coarse_lon[k] = fine_lon[I[k]]
    * J          * idem

    Returns:
    * M3dCoarse * 3d numpy array
    '''
    step = I[3]-I[2]
    jpk, jpj, jpi = MaskCoarse.shape
        
    M3dCoarse=np.zeros((jpk,jpj,jpi),np.float64)
    for jk in range(jpk):
        for ji in range(jpi):
            for jj in range(jpj):
                if MaskCoarse.mask[jk,jj,ji]:
                    i=I[ji]
                    j=J[jj]
                    mask  = Maskfine.mask[jk,j:j+step,i:i+step]
                    values=       M3dfine[jk,j:j+step,i:i+step]
                    area  = Maskfine.area[   j:j+step,i:i+step]
                    M3dCoarse[jk,jj,ji] = weighted_mean(values[mask], area[mask])
    return M3dCoarse


def RSTwriter(outfile, var,rst, TheMask):
    rst[~TheMask.mask] = 1.e+20
    jpk, jpj, jpi = TheMask.shape
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',1)

    TRN   = 'TRN' + var;
    ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = TheMask.ylevels
    ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = TheMask.zlevels
    ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
    ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x') ); ncvar[:] = rst;


    setattr(ncOUT.variables[TRN]   ,'missing_value',1.e+20                             );
    setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
    setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
    ncOUT.close()




I = np.loadtxt('ogstm_South_West_I_indexes.txt', np.int32) # of M_med mesh
J = np.loadtxt('ogstm_South_West_J_indexes.txt', np.int32)

delta = 871 - 722
I = I-delta
I[I<0] = 0



MaskCoarse = Mask("/marconi_work/OGS_dev_0/DEGRADATION_4_70/PREPROC/preproc_meshgen_forcings/mesh_gen/meshmask_470.nc")
Mask16   = Mask("/marconi_scratch/userexternal/gbolzon0/2017_RA16/wrkdir/MODEL/meshmask.nc")
#  correction of spain point, setting tmask=1 where 1/4 mask is 1
orig_tmask=Mask16.mask.copy()
J_corr=[149, 149, 149, 150, 150, 151, 151, 152, 153]
I_corr=[137, 138, 139, 138, 139, 138, 139, 139, 139]

ncorr=len(J_corr)
for i in range(ncorr):  Mask16.mask[:5, J_corr[i], I_corr[i] ] = True


INPUTDIR="/marconi_scratch/userexternal/gbolzon0/REA_KD490/INIT/"
OUTPUTDIR="/marconi_work/OGS_dev_0/DEGRADATION_4_70/PREPROC/preproc_meshgen_forcings/mesh_gen/INIT_INTERP/"

RST_LIST=glob.glob(INPUTDIR +"RST*nc")


for filename in RST_LIST:
    basename=os.path.basename(filename)
    outfile = OUTPUTDIR + basename
    var = basename.rsplit(".")[2]
    print(outfile)
    
    RST = DataExtractor(Mask16, filename, "TRN" + var).values
    values=RST[orig_tmask]
    print("ORIG:         " , values.min(), values.max())
#  correction of spain point, setting values where tmask was 0 #######
    nearest_values= RST[:5,149,140]
    for i in range(ncorr): 
        RST[:5, J_corr[i], I_corr[i] ] = nearest_values
######################################################################

    rst = interp3d(RST, Mask16, MaskCoarse, I, J)
    values = rst[MaskCoarse.mask]
    if np.isnan(values).any():
        raise ValueError
    print("INTERPOLATED: " ,  values.min(), values.max())
    
    
    RSTwriter(outfile, var, rst, MaskCoarse)
