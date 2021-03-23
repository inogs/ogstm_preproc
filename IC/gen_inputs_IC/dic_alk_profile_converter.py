from commons.mask import Mask
from commons.submask import SubMask
from commons.dataextractor import DataExtractor
from basins import V2
import numpy as np
from commons import netcdf4
import netCDF4

#IngvMask = Mask('/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask_CMCCfor_ogstm.nc')
TheMask = Mask('/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask.nc')
filename='/gpfs/scratch/userexternal/gbolzon0/V7C/RHO/rho.2019.nc'
INPUTDIR="/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/IC/inputs/DIC_and_ALK/"
OUTDIR="/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/IC/TRANSITION/inputs/"

jpk, jpj, jpi = TheMask.shape
#JPK, JPJ, JPI = IngvMask.shape
rho =DataExtractor(TheMask,filename, 'rho')


def dumpfile(outfile,M2d,var):
    ncOUT = netCDF4.Dataset(outfile,'w')
    ncOUT.createDimension('nsub',16)
    ncOUT.createDimension('nav_lev',jpk)
    ncvar=ncOUT.createVariable(var,'f',('nsub','nav_lev'))
    ncvar[:]=M2d
    ncOUT.close()
RHO = np.zeros((17,jpk))

for isub, sub in enumerate(V2.Pred):
    print isub
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
