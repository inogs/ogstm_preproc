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
NCin3=NC.Dataset("U20140501-12:00:00.nc","r");
NCin4=NC.Dataset("V20140501-12:00:00.nc","r");
NCin5=NC.Dataset("W20140501-12:00:00.nc","r");
    
#muT               = np.ones((jpk,jpj,jpi),np.double);
#muT[:,:,:]        = (NCin1.variables['muT'][:,:,:]).copy().astype(np.double)
    
e1u               = np.ones((jpj,jpi),np.double);
e1u[:,:]          = (NCin1.variables['e1u'][0,:,:]).copy().astype(np.double)

e2u               = np.ones((jpj,jpi),np.double);
e2u[:,:]          = (NCin1.variables['e2u'][0,:,:]).copy().astype(np.double)

e1v               = np.ones((jpj,jpi),np.double);
e1v[:,:]          = (NCin1.variables['e1v'][0,:,:]).copy().astype(np.double)

e2v               = np.ones((jpj,jpi),np.double);
e2v[:,:]          = (NCin1.variables['e2v'][0,:,:]).copy().astype(np.double)

e1t               = np.ones((jpj,jpi),np.double);
e1t[:,:]          = (NCin1.variables['e1t'][0,:,:]).copy().astype(np.double)

e2t               = np.ones((jpj,jpi),np.double);
e2t[:,:]          = (NCin1.variables['e2t'][0,:,:]).copy().astype(np.double)

# e1f               = np.ones((jpj,jpi),np.double);
# e1f[:,:]          = (NCin1.variables['e1f'][0,:,:]).copy().astype(np.double)

# e2f               = np.ones((jpj,jpi),np.double);
# e2f[:,:]          = (NCin1.variables['e2f'][0,:,:]).copy().astype(np.double)

e3u_0             = np.ones((jpk,jpj,jpi),np.double);
e3u_0[:,:,:]      = (NCin1.variables['e3u_0'][0,:,:,:]).copy().astype(np.double)

e3v_0             = np.ones((jpk,jpj,jpi),np.double);
e3v_0[:,:,:]      = (NCin1.variables['e3v_0'][0,:,:,:]).copy().astype(np.double)

e3w_0             = np.ones((jpk,jpj,jpi),np.double);
e3w_0[:,:,:]      = (NCin1.variables['e3w_0'][0,:,:,:]).copy().astype(np.double)

e3t_0             = np.ones((jpk,jpj,jpi),np.double);
e3t_0[:,:,:]      = (NCin1.variables['e3t_0'][0,:,:,:]).copy().astype(np.double)
    
gdept_0           = np.ones((jpk,jpj,jpi),np.double);
gdept_0[:,:,:]    = (NCin1.variables['gdept_0'][0,:,:,:]).copy().astype(np.double)
    
gdepw_0           = np.ones((jpk,jpj,jpi),np.double);
gdepw_0[:,:,:]    = (NCin1.variables['gdepw_0'][0,:,:,:]).copy().astype(np.double)

mbathy            = np.ones((jpj,jpi),np.int);
mbathy[:,:]       = (NCin1.variables['mbathy'][0,:,:]).copy().astype(np.int)

umask             = np.ones((jpk,jpj,jpi),np.double);
umask[:,:,:]      = (NCin1.variables['umask'][0,:,:,:]).copy().astype(np.double)

vmask             = np.ones((jpk,jpj,jpi),np.double);
vmask[:,:,:]      = (NCin1.variables['vmask'][0,:,:,:]).copy().astype(np.double)

fmask             = np.ones((jpk,jpj,jpi),np.double);
fmask[:,:,:]      = (NCin1.variables['fmask'][0,:,:,:]).copy().astype(np.double)

tmask             = np.ones((jpk,jpj,jpi),np.double);
tmask[:,:,:]      = (NCin1.variables['tmask'][0,:,:,:]).copy().astype(np.double)

ssh               = np.ones((jpj,jpi),np.double);
ssh[:,:]          = (NCin2.variables['sossheig'][0,:,:]).copy().astype(np.double)
    
e3u_t             = np.ones((jpk,jpj,jpi),np.double);
e3u_t[:,:,:]      = (NCin3.variables['e3u'][0,:,:,:]).copy().astype(np.double)
    
e3v_t             = np.ones((jpk,jpj,jpi),np.double);
e3v_t[:,:,:]      = (NCin4.variables['e3v'][0,:,:,:]).copy().astype(np.double)
    
e3w_t             = np.ones((jpk,jpj,jpi),np.double);
e3w_t[:,:,:]      = (NCin5.variables['e3w'][0,:,:,:]).copy().astype(np.double)
    
e3t_t             = np.ones((jpk,jpj,jpi),np.double);
e3t_t[:,:,:]      = (NCin2.variables['e3t'][0,:,:,:]).copy().astype(np.double)

fsdept_n          = np.ones((jpk,jpj,jpi),np.double);
fsdept_n[:,:,:]   = (NCin2.variables['fsdept_n'][0,:,:,:]).copy().astype(np.double)

fsde3w_n          = np.ones((jpk,jpj,jpi),np.double);
fsde3w_n[:,:,:]   = (NCin2.variables['fsde3w_n'][0,:,:,:]).copy().astype(np.double)

fsdepw_n          = np.ones((jpk,jpj,jpi),np.double);
fsdepw_n[:,:,:]   = (NCin5.variables['fsdepw_n'][0,:,:,:]).copy().astype(np.double)


# calculation of new e3 scale factors
# new e3t is calculated following Leclair and Madec (2011) as ratio between e3t_0 and depth H0: 
# e3t is based on e3t* = e3t_0 * (1 + ssh/H0) 
# basically: fse3t_a = e3t_0*(1+ssha/ht_0) as from domvvl.F90
# beware to division by zero: added an IF to calculate only on non-zero depth points
#
# new e3u e3v e3w are computed following the domvvl.F90 (report MERSEA by Levier et al., 2007 is outdated) 

# e3f_new       = np.ones((jpk,jpj,jpi),np.double);
e3t_new       = np.ones((jpk,jpj,jpi),np.double);

for k in range(jpk):
    print 'e3_t for: ', k
    for j in range(jpj):
       for i in range(jpi):
           if e3t_0[0:mbathy[j,i],j,i].sum(0) >0.:
              e3t_new[k,j,i] = e3t_0[k,j,i] * (1.0 +  ssh[j,i]/e3t_0[0:mbathy[j,i],j,i].sum(0) * tmask[k,j,i])
           else:
              e3t_new[k,j,i] = e3t_0[k,j,i]

e3u_new       = np.ones((jpk,jpj,jpi),np.double);
e3v_new       = np.ones((jpk,jpj,jpi),np.double);
e3w_new       = np.ones((jpk,jpj,jpi),np.double);

for k in range(jpk):
    print 'e3_u,_v for: ', k
    for j in range(jpj-1):
        for i in range(jpi-1):
            e3u_new[k,j,i] = 0.5 * (umask[k,j,i]/(e1u[j,i]*e2u[j,i])) * (e1t[j,i]*e2t[j,i]* (e3t_new[k,j,i] - e3t_0[k,j,i]) + e1t[j,i+1]*e2t[j,i+1]*(e3t_new[k,j,i+1]-e3t_0[k,j,i+1]))
            e3v_new[k,j,i] = 0.5 * (vmask[k,j,i]/(e1v[j,i]*e2v[j,i])) * (e1t[j,i]*e2t[j,i]* (e3t_new[k,j,i] - e3t_0[k,j,i]) + e1t[j+1,i]*e2t[j+1,i]*(e3t_new[k,j+1,i]-e3t_0[k,j+1,i]))

for k in range(jpk):
    print 'update e3_u,_v for: ', k
    for j in range(jpj-1):
        for i in range(jpi-1):
           e3u_new[k,j,i] = e3u_0[k,j,i] + e3u_new[k,j,i]
           e3v_new[k,j,i] = e3v_0[k,j,i] + e3v_new[k,j,i]
    e3u_new[k,jpj,:] = e3u_0[k,jpj,:]
    e3v_new[k,jpj,:] = e3v_0[k,jpj,:]
    e3u_new[k,:,jpi] = e3u_0[k,:,jpi]
    e3v_new[k,:,jpi] = e3v_0[k,:,jpi]
          
for j in range(jpj):
    for i in range(jpi):
           e3w_new[0,j,i] = e3w_0[0,j,i] + e3t_new[0,j,i] - e3t_0[0,j,i]

for k in np.arange(1,jpk):
    print 'e3_w for: ', k 
    for j in range(jpj):
       for i in range(jpi):
           e3w_new[k,j,i] = e3w_0[k,j,i] + (1. - 0.5 * tmask[k,j,i]) * (e3t_new[k-1,j,i] - e3t_0[k-1,j,i]) + 0.5 * tmask[k,j,i] * (e3t_new[k,j,i] - e3t_0[k,j,i])

# compute fsdept and fsdepw following dom_vvl_init of domvvl.F90
# this corresponds to computation in report MERSEA

fsdept       = np.ones((jpk,jpj,jpi),np.double);
fsdepw       = np.ones((jpk,jpj,jpi),np.double);
fsde3w       = np.ones((jpk,jpj,jpi),np.double);
zcoef        = np.ones((jpk,jpj,jpi),np.double);

for j in range(jpj):
    for i in range(jpi):
        fsdept[0,j,i] = 0.5 * e3w_new[0,j,i]
        fsdepw[0,j,i] = 0.0
        fsde3w[0,j,i] = fsdept[0,j,i] - ssh[j,i]

for k in np.arange(1,jpk):
    print 'fsdep_t,_w for: ', k 
    for j in range(jpj):
       for i in range(jpi):
# NON ABBIAMO wmask in mesh_mask ...sembra il zcoef serva in punti specifici jpk=mikt ??? ...per ora commento e considero zcoef = 0.
#           zcoef[k,j,i]  = tmask[k,j,i] - wmask[k,j,i]
           fsdepw[k,j,i] = fsdepw[k-1,j,i] + e3t_new[k-1,j,i]
#           fsdept[k,j,i] = zcoef * (fsdepw[k,j,i] + 0.5 * e3w_new[k,j,i]) + (1-zcoef) * (fsdept[k-1,j,i] + e3w_new[k,j,i])
           fsdept[k,j,i] = fsdept[k-1,j,i] + e3w_new[k,j,i]
           fsde3w[k,j,i] = fsdept[k,j,i] - ssh[j,i]



NCin1.close()
NCin2.close()
NCin3.close()
NCin4.close()
NCin5.close()


##############################################################
# write meshmask netcdf file !
##############################################################
ncOUT=NC.Dataset('check_e3all-V2.nc',"w");

ncOUT.createDimension('x',jpi);
ncOUT.createDimension('y',jpj);
ncOUT.createDimension('z',jpk);
    

ncvar    = ncOUT.createVariable('e3t_new' ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3t_new;
ncvar    = ncOUT.createVariable('e3t_t'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3t_t  ;
ncvar    = ncOUT.createVariable('tdiff'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3t_new-e3t_t;

ncvar    = ncOUT.createVariable('e3u_new' ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3u_new;
ncvar    = ncOUT.createVariable('e3u_t'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3u_t  ;
ncvar    = ncOUT.createVariable('udiff'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3u_new-e3u_t;

ncvar    = ncOUT.createVariable('e3v_new' ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3v_new;
ncvar    = ncOUT.createVariable('e3v_t'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3v_t  ;
ncvar    = ncOUT.createVariable('vdiff'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3v_new-e3v_t;

ncvar    = ncOUT.createVariable('e3w_new' ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3w_new;
ncvar    = ncOUT.createVariable('e3w_t'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3w_t  ;
ncvar    = ncOUT.createVariable('wdiff'   ,'d',('z', 'y', 'x'))   ; ncvar[:] = e3w_new-e3w_t;

ncvar    = ncOUT.createVariable('fsdept'  ,'d',('z', 'y', 'x'))   ; ncvar[:] = fsdept;
ncvar    = ncOUT.createVariable('fsdept_n','d',('z', 'y', 'x'))   ; ncvar[:] = fsdept_n;
ncvar    = ncOUT.createVariable('deptdiff','d',('z', 'y', 'x'))   ; ncvar[:] = fsdept-fsdept_n;

ncvar    = ncOUT.createVariable('fsde3w'  ,'d',('z', 'y', 'x'))   ; ncvar[:] = fsde3w;
ncvar    = ncOUT.createVariable('fsde3w_n','d',('z', 'y', 'x'))   ; ncvar[:] = fsde3w_n;
ncvar    = ncOUT.createVariable('de3wdiff','d',('z', 'y', 'x'))   ; ncvar[:] = fsde3w-fsde3w_n;

ncvar    = ncOUT.createVariable('fsdepw'  ,'d',('z', 'y', 'x'))   ; ncvar[:] = fsdepw;
ncvar    = ncOUT.createVariable('fsdepw_n','d',('z', 'y', 'x'))   ; ncvar[:] = fsdepw_n;
ncvar    = ncOUT.createVariable('depwdiff','d',('z', 'y', 'x'))   ; ncvar[:] = fsdepw-fsdepw_n;

ncOUT.close()
    
