# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES


import numpy as np

#import scipy.io.netcdf as NC
import netCDF4 as NC

# Script to create meshmask file

class OrigMask():

    def __init__(self,inputfile):
        NCin=NC.Dataset(inputfile,"r")
        self.Lon = (NCin.variables['XC'][0,:]).copy()
        self.Lat = (NCin.variables['YC'][:,0]).copy()
        self.nav_lev = - (NCin.variables['Z'][:]).copy().astype(np.float);
        self.nc_handler = NCin

    def get_wp(self):
        '''
        Prints some information about the number of waterpooints in longitude
        The aim is to provide some information about where to cut this mesh to get the ogstm meshmask.
        '''
        NCin=self.nc_handler
        jpi=NCin.dimensions['X'].size
        jpj=NCin.dimensions['Y'].size
        jpk=NCin.dimensions['Z'].size
        tmask =np.zeros((jpk,jpj,jpi),dtype=np.int64)
        HFacC = (NCin.variables['HFacC'][:,:,:]).copy() #shape is (121, 380, 1307)
        idx_wet = HFacC > 0.
        tmask[idx_wet] = 1

        waterpoints_longitude = tmask.sum(axis=0).sum(axis=0)

        Gibraltar_lon = -5.75
        imed = self.Lon>Gibraltar_lon
        med_waterpoints = waterpoints_longitude[imed].sum()
        print "mediterranean waterpoints:", med_waterpoints
        print "atlantic waterpoints:"
        for lon in np.arange(-9,-8,1./16):
            iatl = (self.Lon < Gibraltar_lon) & (self.Lon> lon)
            atl_waterpoints = waterpoints_longitude[iatl].sum()
            print lon, atl_waterpoints

    def getCutLocation(self,lon):
        ind =np.argmin(np.abs(self.Lon-lon))
        return ind



time = 1;
z_a  = 1;
y_a  = 1;
x_a  = 1;




def create_meshmask_nc(OrigMaskobj,outfile,free_surface=False):
    
    NCin=OrigMaskobj.nc_handler

    jpi=NCin.dimensions['X'].size
    jpj=NCin.dimensions['Y'].size
    jpk=NCin.dimensions['Z'].size
    time = 1

#    double fmask(time, z, y, x) ;
    fmask = np.zeros((time,jpk,jpj,jpi),np.double);
    print("fmask not used, assigned to 0.")
#   fmask[0,:,:,:]=(NCin.variables['fmask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)

#    fmask[0,:,1, 0]    = 0.;
#    fmask[0,:,0, 1]    = 0.;
#    fmask[0,:,-1,0]    = 0.;
#    fmask[0,:,0,-1]    = 0.;

#    double gdept(time, z, y_a, x_a) ;
    gdept              = np.zeros((time,jpk,y_a,x_a),np.double);
    gdept[0,0:jpk,0,0] = - (NCin.variables['Z'][:jpk]).copy().astype(np.double);

#    double gdepw(time, z, y_a, x_a) ;
    gdepw              = np.zeros((time,jpk,y_a,x_a),np.double);
    gdepw[0,0:jpk,0,0] = - (NCin.variables['Zp1'][0:jpk]).copy().astype(np.double);
    

#    double glamt(time, z_a, y, x) ;
    glamt = np.ones((time,z_a,jpj,jpi),np.double);
    glamt[0,0,:,:]=(NCin.variables['XC'][:,:]).copy().astype(np.double);
        
#    double glamu(time, z_a, y, x) ;
    glamu = np.ones((time,z_a,jpj,jpi),np.double);
    glamu[0,0,:,:]=(NCin.variables['XG'][0:jpj,1:jpi+1]).copy().astype(np.double);
        
#    double glamv(time, z_a, y, x) ;
    glamv = np.ones((time,z_a,jpj,jpi),np.double);
    glamv[0,0,:,:]=(NCin.variables['XC'][:,:]).copy().astype(np.double);

#    double glamf(time, z_a, y, x) ;
    glamf = np.ones((time,z_a,jpj,jpi),np.double);
    glamf[0,0,:,:]=(NCin.variables['XG'][1:jpj+1,1:jpi+1]).copy().astype(np.double);

        
#    double gphit(time, z_a, y, x) ;
    gphit = np.ones((time,z_a,jpj,jpi),np.double);
    gphit[0,0,:,:]=(NCin.variables['YC'][:,:]).copy().astype(np.double);
        
#    double gphiu(time, z_a, y, x) ;
    gphiu = np.ones((time,z_a,jpj,jpi),np.double);
    gphiu[0,0,:,:]=(NCin.variables['YC'][:,:]).copy().astype(np.double);
    
#    double gphiv(time, z_a, y, x) ;
    gphiv = np.ones((time,z_a,jpj,jpi),np.double);
    gphiv[0,0,:,:]=(NCin.variables['YG'][1:jpj+1,0:jpi]).copy().astype(np.double);
        
#    double gphif(time, z_a, y, x) ;
    gphif = np.ones((time,z_a,jpj,jpi),np.double);
    gphif[0,0,:,:]=(NCin.variables['YG'][1:jpj+1,1:jpi+1]).copy().astype(np.double);
    
#    double ff(time, z_a, y, x) ;
    ff = np.zeros((time,z_a,jpj,jpi),np.double);
#   ff[0,0,:,:]=(NCin.variables['ff'][0,:,lon_cut:]).copy().astype(np.double);


#    float nav_lat(y, x) ;
    nav_lat = np.ones((jpj,jpi),np.float);
    nav_lat[:,:] = (NCin.variables['YC'][:,:]).copy().astype(np.float);
    
#    float nav_lev(z) ;
    nav_lev = np.ones((jpk,),np.float);
    nav_lev[:] = -(NCin.variables['Z'][:]).copy().astype(np.float);

#    float nav_lon(y, x) ;
    nav_lon = np.ones((jpj,jpi),np.float);
    nav_lon[:,:] = (NCin.variables['XC'][:,:]).copy().astype(np.float);
    
    
#    double e1f(time, z_a, y, x) ;
    e1f = np.ones((time,z_a,jpj,jpi),np.double);
    e1f[0,0,:,:] =  (NCin.variables['dxV'][1:jpj+1,1:jpi+1]).copy().astype(np.double)

#    double e1t(time, z_a, y, x) ;
    e1t = np.ones((time,z_a,jpj,jpi),np.double);
    e1t[0,0,:,:] =  (NCin.variables['dxF'][:,:]).copy().astype(np.double)
    
#    double e1u(time, z_a, y, x) ;
    e1u = np.ones((time,z_a,jpj,jpi),np.double);
    e1u[0,0,:,:] =  (NCin.variables['dxC'][:,1:jpi+1]).copy().astype(np.double)

#    double e1v(time, z_a, y, x) ;
    e1v = np.ones((time,z_a,jpj,jpi),np.double);
    e1v[0,0,:,:] =  (NCin.variables['dxG'][1:jpj+1,:]).copy().astype(np.double)
    
    
#    double e2f(time, z_a, y, x) ;
    e2f = np.ones((time,z_a,jpj,jpi),np.double);
    e2f[0,0,:,:] =  (NCin.variables['dyU'][1:jpj+1,1:jpi+1]).copy().astype(np.double)
    
    
#    double e2t(time, z_a, y, x) ;
    e2t = np.ones((time,z_a,jpj,jpi),np.double)
    e2t[0,0,:,:] =  (NCin.variables['dyF'][:,:]).copy().astype(np.double)
    
#    double e2u(time, z_a, y, x) ;
    e2u = np.ones((time,z_a,jpj,jpi),np.double);
    e2u[0,0,:,:] =  (NCin.variables['dyG'][:,1:jpi+1]).copy().astype(np.double)
    
#    double e2v(time, z_a, y, x) ;
    e2v = np.ones((time,z_a,jpj,jpi),np.double);
    e2v[0,0,:,:] =  (NCin.variables['dyC'][1:jpj+1,:]).copy().astype(np.double)
    
    e3t_1d=np.zeros((jpk))
    Zp1 = - (NCin.variables['Zp1'][:]).copy().astype(np.double)
    e3t_1d=np.diff(Zp1)
    
    e3w_1d=np.zeros((jpk))
    Z = - (NCin.variables['Z'][:]).copy().astype(np.double)
    
    e3w_1d[0]=Z[0]*2.
    for k in range(1,jpk):
        e3w_1d[k] = Z[k] -Z[k-1]


    if free_surface:

#    double e3t(time, z, y, x) ;
        e3t = np.ones((time,jpk,jpj,jpi),np.double)
        HFacC =  (NCin.variables['HFacC'][:,:,:]).copy().astype(np.double)
        for jj in range(jpj):
            for ji in range(jpi):
                e3t[0,:,jj,ji] =  e3t_1d[:] * HFacC[:,jj,ji]

#    double e3u(time, z, y, x) ;
        e3u = np.ones((time,jpk,jpj,jpi),np.double);
        HFacW =  (NCin.variables['HFacW'][:,:,1:jpi+1]).copy().astype(np.double)
        for jj in range(jpj):
            for ji in range(jpi):
                e3u[0,:,jj,ji] =  e3t_1d[:] * HFacW[:,jj,ji]

#    double e3v(time, z, y, x) ;
        e3v = np.ones((time,jpk,jpj,jpi),np.double);
        HFacS =  (NCin.variables['HFacS'][:,1:jpj+1,:]).copy().astype(np.double)
        for jj in range(jpj):
            for ji in range(jpi):
                e3v[0,:,jj,ji] =  e3t_1d[:] * HFacS[:,jj,ji]

#    double e3w(time, z, y, x) ;
        e3w = np.ones((time,jpk,jpj,jpi),np.double);
        for jj in range(jpj):
            for ji in range(jpi):
                e3w[0,:,jj,ji] =  e3w_1d[:] * HFacC[:,jj,ji]
    else:
#    double e3t(time, z, y_a, x_a) ;
        e3t            = np.ones((time,jpk,y_a,x_a),np.double);
        e3t[0,:,0,0]   =  e3t_1d[:]
    
#    double e3w(time, z, y_a, x_a) ;
        e3w            = np.zeros((time,jpk,y_a,x_a),np.double);
        e3w[0,:,0,0]   = e3w_1d[:]

#    double tmask(time, z, y, x) ;

    tmask =np.zeros((time,jpk,jpj,jpi),dtype=np.double)
    HFacC = (NCin.variables['HFacC'][:,:,:]).copy() #shape is (121, 380, 1307)
    idx_wet = HFacC > 0.
    tmask[0,idx_wet] = 1.

    tmask[0,:,0, :] =0.;
    tmask[0,:,:, 0] =0.;
    tmask[0,:,jpj-1,:] =0.;
    tmask[0,:,:,jpi-1] =0.;
    tmask[0,jpk-1,:,:] =0.;


#    double umask(time, z, y, x) ;
    umask = np.zeros((time,jpk,jpj,jpi),np.double);
    HFacW = (NCin.variables['HFacW'][:,:,1:jpi+1]).copy() #shape is (121, 380, 1307)
    idx_wet = HFacW > 0.
    umask[0,idx_wet] = 1.

    umask[0,:,0, :] =0.;
    umask[0,:,:, 0] =0.;
    umask[0,:,jpj-1,:] =0.;
    umask[0,:,:,jpi-1] =0.;
#   umask[0,:,:,-2] =0.; To be checked !!
    umask[0,jpk-1,:,:] =0.;

#    double vmask(time, z, y, x) ;

    vmask = np.ones((time,jpk,jpj,jpi),np.double);
    HFacS = (NCin.variables['HFacS'][:,1:jpj+1,:]).copy() #shape is (121, 380, 1307)
    idx_wet = HFacS > 0.
    vmask[0,idx_wet] = 1.
    vmask[0,:,:, 0] =0.;
    vmask[0,:,0, :] =0.;
    vmask[0,:,jpj-1,:] =0.;
    vmask[0,:,:,jpi-1] =0.;
#    vmask[0,:,-2,:] =0.; To be checked !!
    vmask[0,jpk-1,:,:] =0.;

    ##############################################################
    # write meshmask netcdf file !
    ##############################################################
    ncOUT=NC.Dataset(outfile,"w");

    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',time)
    
    ncOUT.createDimension('x_a',x_a);
    ncOUT.createDimension('y_a',y_a);
    ncOUT.createDimension('z_a',z_a);

#   ncvar    = ncOUT.createVariable('coastp','d',('y','x'))                  ; ncvar[:] = coastp;  
    ncvar    = ncOUT.createVariable('e1f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e1f   ;
    ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e1t   ;
    ncvar    = ncOUT.createVariable('e1u'   ,'d',('time','z_a', 'y', 'x')   ); ncvar[:] = e1u   ;
    ncvar    = ncOUT.createVariable('e1v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e1v   ;
    ncvar    = ncOUT.createVariable('e2f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e2f   ;
    ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e2t   ;
    ncvar    = ncOUT.createVariable('e2u'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2u   ;
    ncvar    = ncOUT.createVariable('e2v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2v   ;     
    if free_surface:   
        ncvar    = ncOUT.createVariable('e3t_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3t   ;
        ncvar    = ncOUT.createVariable('e3u_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3u   ;
        ncvar    = ncOUT.createVariable('e3v_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3v   ;
        ncvar    = ncOUT.createVariable('e3w_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3w   ;   
    else:
        ncvar    = ncOUT.createVariable('e3t'   ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[:] = e3t   ;
        ncvar    = ncOUT.createVariable('e3w'   ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[:] = e3w   ;          
    ncvar    = ncOUT.createVariable('ff'    ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = ff    ;      
    ncvar    = ncOUT.createVariable('fmask' ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = fmask ;    
    ncvar    = ncOUT.createVariable('gdept' ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[:] = gdept ;
    ncvar    = ncOUT.createVariable('gdepw' ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[:] = gdepw ;
    ncvar    = ncOUT.createVariable('glamf'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamf ;     
    ncvar    = ncOUT.createVariable('glamt'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamt ;
    ncvar    = ncOUT.createVariable('glamu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamu ;     
    ncvar    = ncOUT.createVariable('glamv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamv ;
    ncvar    = ncOUT.createVariable('gphif'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphif ;     
    ncvar    = ncOUT.createVariable('gphit'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphit ;
    ncvar    = ncOUT.createVariable('gphiu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphiu ;     
    ncvar    = ncOUT.createVariable('gphiv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphiv ;
    ncvar    = ncOUT.createVariable('nav_lat','f',('y','x'))                 ; ncvar[:] = nav_lat;
    ncvar    = ncOUT.createVariable('nav_lev' ,'f',('z',))                   ; ncvar[:] = nav_lev;
    ncvar    = ncOUT.createVariable('nav_lon','f',('y','x'))                 ; ncvar[:] = nav_lon;
#	float time(time) ;
#	short time_steps(time) ;
    ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = tmask 
    ncvar    = ncOUT.createVariable('umask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = umask 
    ncvar    = ncOUT.createVariable('vmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = vmask 
    ncOUT.close()
    
