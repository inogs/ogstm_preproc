# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES


import numpy as np

#import scipy.io.netcdf as NC
import netCDF4 as NC

# Script to create meshmask file

    
#NCin1=NC.Dataset("prova_mesh.nc","r");
NCin1=NC.Dataset("mesh_mask_cut.nc","r");

jpi=NCin1.dimensions['x'].size
jpj=NCin1.dimensions['y'].size
jpk=NCin1.dimensions['z'].size

#NCin2=NC.Dataset("/gpfs/scratch/userexternal/plazzari/eas_v12/FORCINGS/OUTPUT/T20140101-12:00:00.nc","r");
#NCin2=NC.Dataset("T20140101-12:00:00.nc","r");
#NCin2=NC.Dataset("T20140201-12:00:00.nc","r");
NCin2=NC.Dataset("T20140501-12:00:00.nc","r");
#NCin2=NC.Dataset("T20141001-12:00:00.nc","r");
    
#muT               = np.ones((jpk,jpj,jpi),np.double);
#muT[:,:,:]        = (NCin1.variables['muT'][:,:,:]).copy().astype(np.double)
    
e3t_0             = np.ones((jpk,jpj,jpi),np.double);
e3t_0[:,:,:]      = (NCin1.variables['e3t_0'][0,:,:,:]).copy().astype(np.double)
    
gdept_0           = np.ones((jpk,jpj,jpi),np.double);
gdept_0[:,:,:]    = (NCin1.variables['gdept_0'][0,:,:,:]).copy().astype(np.double)

mbathy            = np.ones((jpj,jpi),np.int);
mbathy[:,:]       = (NCin1.variables['mbathy'][0,:,:]).copy().astype(np.int)

tmask             = np.ones((jpk,jpj,jpi),np.double);
tmask[:,:,:]      = (NCin1.variables['tmask'][0,:,:,:]).copy().astype(np.double)

ssh               = np.ones((jpj,jpi),np.double);
ssh[:,:]          = (NCin2.variables['sossheig'][0,:,:]).copy().astype(np.double)
    
e3t_t             = np.ones((jpk,jpj,jpi),np.double);
e3t_t[:,:,:]      = (NCin2.variables['e3t'][0,:,:,:]).copy().astype(np.double)

e3t_new       = np.ones((jpk,jpj,jpi),np.double);
for k in range(jpk):
    print k
    for j in range(jpj):
       for i in range(jpi):
           if e3t_0[0:mbathy[j,i],j,i].sum(0) >0:
              e3t_new[k,j,i] = e3t_0[k,j,i] * (1 + ssh[j,i]/e3t_0[0:mbathy[j,i],j,i].sum(0) * tmask[k,j,i])
           else:
              e3t_new[k,j,i] = e3t_t[k,j,i]
          
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
    
