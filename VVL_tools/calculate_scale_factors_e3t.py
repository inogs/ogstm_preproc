# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES


import numpy as np

#import scipy.io.netcdf as NC
import netCDF4 as NC

# Script to create meshmask file

    
NCin1=NC.Dataset("prova_mesh.nc","r");

jpi=NCin1.dimensions['x'].size
jpj=NCin1.dimensions['y'].size
jpk=NCin1.dimensions['z'].size

NCin2=NC.Dataset("/gpfs/scratch/userexternal/plazzari/eas_v12/FORCINGS/OUTPUT/T20140101-12:00:00.nc","r");
    
muT               = np.ones((jpk,jpj,jpi),np.double);
muT[:,:,:]        = (NCin1.variables['muT'][:,:,:]).copy().astype(np.double)
    
e3t               = np.ones((jpk,jpj,jpi),np.double);
e3t[:,:,:]        = (NCin1.variables['e3t'][0,:,:,:]).copy().astype(np.double)
    
ssh               = np.ones((jpj,jpi),np.double);
ssh[:,:]          = (NCin2.variables['sossheig'][0,:,:]).copy().astype(np.double)
    
e3t_t             = np.ones((jpk,jpj,jpi),np.double);
e3t_t[0:125,:,:]  = (NCin2.variables['e3t'][0,:,:,:]).copy().astype(np.double)


e3t_new       = np.ones((jpk,jpj,jpi),np.double);
for k in range(jpk):
    for j in range(jpj):
       for i in range(jpi):
           e3t_new[k,j,i] = e3t[k,j,i] + ssh[j,i]*muT[k,j,i]
          
NCin1.close()
NCin2.close()


##############################################################
# write meshmask netcdf file !
##############################################################
ncOUT=NC.Dataset('check.nc',"w");

ncOUT.createDimension('x',jpi);
ncOUT.createDimension('y',jpj);
ncOUT.createDimension('z',jpk);
    

ncvar    = ncOUT.createVariable('e3t_new','d',('z', 'y', 'x'))   ; ncvar[:] = e3t_new;
ncvar    = ncOUT.createVariable('e3t_t'  ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3t_t  ;
ncvar    = ncOUT.createVariable('diff'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3t_new-e3t_t;
ncOUT.close()
    
