import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Intepolates climatology from 1 km to spatial resolution defined by meshmask.nc
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--inputfile', '-i',
                                type = str,
                                required = True,
                                help = '''Clim file at 1 km'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '''maskfile of the output climatology'''
                                )
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = True,
                                help = ''' path of the output climatology'''
                                )
    return parser.parse_args()

args = argument()

from commons.mask import Mask
from Sat import interp2d 
from Sat import SatManager as Sat
import numpy as np
import netCDF4

TheMask=Mask(args.maskfile)
CLIM_FILE=args.inputfile
MEAN,STD = Sat.readClimatology(CLIM_FILE)

jpk,jpj,jpi = TheMask.shape
x = TheMask.xlevels[0,:]
y = TheMask.ylevels[:,0]

x1km = Sat.masks.KD490mesh.lon
y1km = Sat.masks.KD490mesh.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

CLIM = np.zeros((365, jpj,jpi), dtype=[('MEAN',np.float32)])

for julian in range(365):
    print julian
    Mfine=MEAN[julian,:,:]
    M,NP  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END)
    CLIM[julian,:,:] = M

ncOUT = netCDF4.Dataset(args.outfile,'w')
ncOUT.createDimension('lon',jpi)
ncOUT.createDimension('lat',jpj)
ncOUT.createDimension('time',365)
ncvar=ncOUT.createVariable('Mean','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIM['MEAN']
ncOUT.close()


