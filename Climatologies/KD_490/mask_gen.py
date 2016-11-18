'''
Saves the mask file for KD490 1km files
Executed once and for all.
'''

import numpy as np
import scipy.io.netcdf as NC
COUNT = np.load('Kd_map_occurency.npy')


filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/20150831_d-OC_CNR-L3-KD490-MedOC4AD4_SAM_1KM-MED-REP-v01.nc"

ncIN = NC.netcdf_file(filename,'r')
LON=ncIN.variables['lon'].data.copy()
LAT=ncIN.variables['lat'].data.copy()
ncIN.close()

tmask = COUNT > 0
jpj,jpi=COUNT.shape
ncOUT = NC.netcdf_file('KD490_1km_meshmask.nc','w')
ncOUT.createDimension("lon", jpi)
ncOUT.createDimension("lat", jpj)
ncvar=ncOUT.createVariable("lon", "f", ('lon',))
ncvar[:]=LON
ncvar=ncOUT.createVariable("lat", "f", ('lat',))
ncvar[:]=LAT
ncvar=ncOUT.createVariable('tmask', 'b', ('lat','lon'))
ncvar[:]=tmask
ncOUT.close()


ncOUT = NC.netcdf_file('KD490_occurrency.nc','w')
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




