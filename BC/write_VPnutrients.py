import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates VPnutrients.nc, with a single season
    starting from profiles in inputdir, in NetCDF files called
    N1p.nc, N3n.nc etc.
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'Input directory with profile files for each variable'
                                )
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )
    return parser.parse_args()

args = argument()

from commons.mask import Mask
from commons import netcdf4
import netCDF4
from commons.utils import addsep
TheMask=Mask(args.maskfile, loadtmask=False)
INPUTDIR=addsep(args.inputdir)

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


