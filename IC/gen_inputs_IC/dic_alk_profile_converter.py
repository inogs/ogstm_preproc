from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons import density
from bitsea.basins import V2
import numpy as np
from bitsea.commons import netcdf4
import netCDF4

MaskCMCC = Mask("/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask_CMCC.nc")
TheMask = Mask("/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc")
filename='/g100_work/OGS_devC/V11C/ICs/RUN_QUID/IC/T.transition.nc'
INPUTDIR="/g100_work/OGS_devC/V11C/ICs/RUN_QUID/IC/inputs/"
OUTDIR=INPUTDIR

jpk, jpj, jpi = TheMask.shape
JPK, JPJ, JPI = MaskCMCC.shape
rho = density.get_density(filename, MaskCMCC)
rho=rho[:jpk,:,JPI-jpi:]


def dumpfile(outfile,M2d,var):
    ncOUT = netCDF4.Dataset(outfile,'w')
    ncOUT.createDimension('nsub',16)
    ncOUT.createDimension('nav_lev',jpk)
    ncvar=ncOUT.createVariable(var,'f',('nsub','nav_lev'))
    ncvar[:]=M2d
    ncOUT.close()
RHO = np.zeros((17,jpk))

for isub, sub in enumerate(V2.Pred):
    print(isub)
    S = SubMask(sub,maskobject=TheMask)    
    for k in range(jpk):
        submask = S.mask[k,:,:]
        V=rho[k,:,:]
        n = submask.sum()
        if n> 0:
            Conc = V[submask]
            Weight=TheMask.area[submask]
            Weight_sum      = Weight.sum()
            Mass            = (Conc * Weight).sum()
            Weighted_Mean   = Mass/Weight_sum
            RHO[isub,k] = Weighted_Mean
        else:
            RHO[isub,k] = Weighted_Mean

DIC=netcdf4.readfile(INPUTDIR + "DIC.nc",'DIC')
ALK=netcdf4.readfile(INPUTDIR + "ALK.nc",'ALK')
O3c = RHO[:16,:] * DIC * 12 /1000
O3h = RHO[:16,:] * ALK /1000


dumpfile(OUTDIR+'O3c.nc', O3c, 'O3c')
dumpfile(OUTDIR+'O3h.nc', O3h, 'O3h')
