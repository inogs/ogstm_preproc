import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Step 2 for the generation of a climatology

    Generates the mask file for 1km Sat dataset (CHL or KD490)
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--inputfile', '-i',
                                type = str,
                                required = True,
                                help = '''map_occurency.npy'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '''CHL_1km_meshmask.nc, output file'''
                                )
    parser.add_argument(   '--occurrency', '-o',
                                type = str,
                                required = True,
                                help = '''CHL_occurrency.nc'''
                                )
    return parser.parse_args()

args = argument()



import numpy as np
from netCDF4 import Dataset
COUNT = np.load(args.inputfile)


filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG/20061115_d-OC_CNR-L3-CHL-MedOC4AD4_SAM_1KM-MED-REP-v02.nc"

ncIN = Dataset(filename, 'r')
LON=ncIN.variables['lon'][:]
LAT=ncIN.variables['lat'][:]
ncIN.close()

tmask = COUNT > 0
jpj,jpi=COUNT.shape
ncOUT = Dataset(args.maskfile,'w')
ncOUT.createDimension("lon", jpi)
ncOUT.createDimension("lat", jpj)
ncvar=ncOUT.createVariable("lon", "f", ('lon',))
ncvar[:]=LON
ncvar=ncOUT.createVariable("lat", "f", ('lat',))
ncvar[:]=LAT
ncvar=ncOUT.createVariable('tmask', 'b', ('lat','lon'))
ncvar[:]=tmask
ncOUT.close()


ncOUT = Dataset(args.occurrency,'w')
ncOUT.createDimension("lon", jpi)
ncOUT.createDimension("lat", jpj)
ncvar=ncOUT.createVariable("lon", "f", ('lon',))
ncvar[:]=LON
ncvar=ncOUT.createVariable("lat", "f", ('lat',))
ncvar[:]=LAT
ncvar=ncOUT.createVariable('nData', 'i', ('lat','lon'))
setattr(ncvar,'missing_value',0)
ncvar[:]=COUNT
ncOUT.close()




