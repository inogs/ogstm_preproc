import numpy as np
import netCDF4
from commons.mask import Mask
from commons import netcdf4
from basins import V2

TheMask=Mask('/gss/gss_work/DRES_OGS_BiGe/ateruzzi/RA_24/input/setup/PREPROC/MASK/gdept_3d/ogstm/meshmask.nc',loadtmask=False)
jpk,_,_ = TheMask.shape
FLOAT_CLIM = netcdf4.readfile("/g100/home/userexternal/camadio0/float_preproc/CLIMATOLOGIE/FLOAT_CLIMATOLOGIES/new_clim_SF_V8C_2012_2021_swm2/yr_Avg_O2o.nc", "O2o")
EMODNET    = netcdf4.readfile("/gss/gss_work/DRES_OGS_BiGe/ateruzzi/RA_24/input/setup/PREPROC/IC/inputs/O2o.nc","O2o")

def dumpfile(outfile,M2d,var):
    ncOUT = netCDF4.Dataset(outfile,'w')
    ncOUT.createDimension('nsub',16)
    ncOUT.createDimension('nav_lev',jpk)
    ncvar=ncOUT.createVariable(var,'f',('nsub','nav_lev'))
    ncvar[:]=M2d
    ncOUT.close()


z_layer=np.arange(0,200,10)
z_layer = np.append(z_layer , np.arange(200, 600, 40))
z_layer = np.append(z_layer , np.arange(600,1001, 50))

z_interp = (z_layer[:-1]+z_layer[1:])/2 


nsub,nlev = FLOAT_CLIM.shape
OUT = np.zeros((nsub,jpk),np.float32)

jk_float_bottom = TheMask.getDepthIndex(z_interp[-1])
deep = TheMask.zlevels > 1000

for isub, sub in enumerate(V2.Pred):
    if sub.name=='atl' : continue
    if sub in [V2.aeg, V2.adr1, V2.adr2]:
        OUT[isub,:] = EMODNET[isub,:]
    else:
        profile = FLOAT_CLIM[isub,:]
        nans = np.isnan(profile)
        shift =  FLOAT_CLIM[isub,-1] - EMODNET[isub,jk_float_bottom] 
        print(shift)
        interpolated = np.interp(TheMask.zlevels,z_interp[~nans],profile[~nans])
        interpolated[deep] = EMODNET[isub, deep] + shift
        
        OUT[isub,:] = interpolated

dumpfile("O2o.nc", OUT, "O2o")

