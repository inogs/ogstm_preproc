from commons.mask import Mask
from commons import netcdf4
import netCDF4
TheMask=Mask('/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc', loadtmask=False)
INPUTDIR="/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/inputs/20201219/"
jpk,_,_ = TheMask.shape

ncOUT = netCDF4.Dataset('VPnutrients.nc','w')
ncOUT.createDimension('lev1',jpk)
ncOUT.createDimension('lev2',jpk)
ncOUT.createDimension('seas',1)

ncvar = ncOUT.createVariable('lev1','f',('lev1',))
ncvar[:] = TheMask.zlevels

ncvar = ncOUT.createVariable('lev2','f',('lev2',))
ncvar[:] = TheMask.zlevels

for var in ['N1p','N3n','N5s','O2o']:
    ncvar = ncOUT.createVariable(var,'f',('seas','lev1',))
    ncvar[:]=netcdf4.readfile(INPUTDIR + var + ".nc",var)

ncvar = ncOUT.createVariable('ALK','f',('lev2',))
ncvar[:] = netcdf4.readfile(INPUTDIR + "O3h.nc",'O3h')

ncvar = ncOUT.createVariable('DIC','f',('lev2',))
ncvar[:] = netcdf4.readfile(INPUTDIR + "O3c.nc",'O3c')

ncOUT.close()


