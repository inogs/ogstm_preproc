# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES


import numpy as np
from commons.mask import Mask
import netCDF4 as NC

maskfile  = 'meshmask.nc'

outfile = 'RESTARTS/RST.19950101-00:00:00.O5c.nc'

TheMask=Mask(maskfile)


jpk, jpj, jpi=TheMask.shape

time=1

nav_lev=TheMask.zlevels

#    double PIC(time, z, y, x) ;
PIC = np.ones((time,jpk,jpj,jpi),np.double)*1.e+20;
for jk in range(jpk):
    for jj in range(jpj):
        for ji in range(jpi):
            if TheMask.mask[jk,jj,ji]:
                if nav_lev[jk] < 150 :
                    PIC[0,jk,jj,ji]  = 2.0
                elif nav_lev[jk] >= 150 and nav_lev[jk] < 300:
                    PIC[0,jk,jj,ji]  = 1.0
                else:
                    PIC[0,jk,jj,ji]  = 1.e-11

##############################################################
# write meshmask netcdf file !
##############################################################
var="O5c"
rst=PIC
ncOUT=NC.Dataset(outfile,"w");
ncOUT.createDimension('x',jpi);
ncOUT.createDimension('y',jpj);
ncOUT.createDimension('z',jpk);
ncOUT.createDimension('time',time)
    


TRB   = 'TRB' + var;
TRN   = 'TRN' + var;
ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = TheMask.xlevels
ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = TheMask.ylevels
ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = nav_lev;
ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
ncvar = ncOUT.createVariable(TRB       ,'d',('time','z','y','x')); ncvar[:] = rst;   
ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x')); ncvar[:] = rst; 

setattr(ncOUT.variables[TRB]   ,'missing_value',1e+20                              );     
setattr(ncOUT.variables[TRN]   ,'missing_value',1e+20                              );
setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
ncOUT.close()

