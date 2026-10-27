import numpy as np
import netCDF4 as NC

# Script to create meshmask file
# this version uses two files, mask_file and domain_file,
# because, CMCC's fault, mask file has all tmask, umask, etc.
# domain_file has no tmask etc. but has glamt, gphit, etc.
# with no zeros on land only domains from the domain decomposition

class OrigMask():

    def __init__(self, mask_file, domain_file):
        NCin_msk=NC.Dataset(mask_file,"r")
        NCin_dmn=NC.Dataset(domain_file,"r")
        self.Lon = (NCin_dmn.variables['glamt'][0,:,:]).copy()
        self.Lat = (NCin_dmn.variables['gphit'][0,:,:]).copy()
        self.nav_lev = (NCin_dmn.variables['nav_lev'][:]).copy().astype(float);
        self.nc_handler_msk = NCin_msk
        self.nc_handler_dmn = NCin_dmn


    def get_wp(self):
        '''
        Prints some information about the number of waterpooints in longitude
        The aim is to provide some information about where to cut this mesh to get the ogstm meshmask.
        '''
        NCin_msk=self.nc_handler_msk
        jpi=NCin_msk.dimensions['x'].size
        jpj=NCin_msk.dimensions['y'].size
        tmask= (NCin_msk.variables['tmask'][0,:,:,:]).copy().astype(np.int64) #shape is (141, 715, 1275)
        # Biscay correction
        for j in range(jpj):
            for i in range(jpi):
                if(self.Lon[j,i]<0.0 and self.Lat[j,i]>42):
                    tmask[:,j,i]=0.
        waterpoints_longitude = tmask.sum(axis=0).sum(axis=0)

        Gibraltar_lon = -5.75
        imed = self.Lon>Gibraltar_lon
        med_waterpoints = waterpoints_longitude[imed].sum()
        print("mediterranean waterpoints:", med_waterpoints)
        print("atlantic waterpoints:")
        for lon in np.arange(-9,-8,1./16):
            iatl = (self.Lon < Gibraltar_lon) & (self.Lon> lon)
            atl_waterpoints = waterpoints_longitude[iatl].sum()
            print(lon, atl_waterpoints)

    def getCutLocation(self, lon, lat):
        #ind =np.argmin(np.abs(self.Lon-lon))
        #print(f"cut_location: {ind}")
        #return ind
        dist = (self.Lon - lon)**2 + (self.Lat - lat)**2
        j, i = np.unravel_index(np.nanargmin(dist), dist.shape)
        return j, i



time = 1;
z_a  = 1;
y_a  = 1;
x_a  = 1;




def create_meshmask_nc(OrigMaskobj, outfile, lon_cut, depth_cut, biscay_land=True):
    
    NCin_msk=OrigMaskobj.nc_handler_msk
    NCin_dmn=OrigMaskobj.nc_handler_dmn

    jpi=NCin_dmn.dimensions['x'].size-lon_cut
    jpj=NCin_dmn.dimensions['y'].size
    jpk=NCin_dmn.dimensions['nav_lev'].size-depth_cut #116

#    double coastp(y, x) ;
    coastp = np.zeros((jpj,jpi),np.double)
#    double fmask(time, z, y, x) ;
    fmask = np.ones((time,jpk,jpj,jpi),np.double);
    fmask[0,:,:,:]=(NCin_msk.variables['fmask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)

#    fmask[0,:,1, 0]    = 0.;
#    fmask[0,:,0, 1]    = 0.;
#    fmask[0,:,-1,0]    = 0.;
#    fmask[0,:,0,-1]    = 0.;

#    double gdept(time, z, y, x) ;
    gdept              = np.zeros((time,jpk,jpj,jpi),np.double);
    gdept[0,0:jpk,:,:] = (NCin_msk.variables['gdept_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double);

#    double gdepw(time, z, y_a, x_a) ;
    gdepw              = np.zeros((time,jpk,y_a,x_a),np.double);
    gdepw[0,0:jpk,0,0] = (NCin_msk.variables['gdepw_1d'][0,:jpk]).copy().astype(np.double);
    

#    double glamt(time, z_a, y, x) ;
    glamt = np.ones((time,z_a,jpj,jpi),np.double);
    glamt[0,0,:,:]=(NCin_dmn.variables['glamt'][0,:,lon_cut:]).copy().astype(np.double);
        
#    double glamu(time, z_a, y, x) ;
    glamu = np.ones((time,z_a,jpj,jpi),np.double);
    glamu[0,0,:,:]=(NCin_dmn.variables['glamu'][0,:,lon_cut:]).copy().astype(np.double);
        
#    double glamv(time, z_a, y, x) ;
    glamv = np.ones((time,z_a,jpj,jpi),np.double);
    glamv[0,0,:,:]=(NCin_dmn.variables['glamv'][0,:,lon_cut:]).copy().astype(np.double);

#    double glamf(time, z_a, y, x) ;
    glamf = np.ones((time,z_a,jpj,jpi),np.double);
    glamf[0,0,:,:]=(NCin_dmn.variables['glamf'][0,:,lon_cut:]).copy().astype(np.double);
        
#    double gphit(time, z_a, y, x) ;
    gphit = np.ones((time,z_a,jpj,jpi),np.double);
    gphit[0,0,:,:]=(NCin_dmn.variables['gphit'][0,:,lon_cut:]).copy().astype(np.double);
        
#    double gphiu(time, z_a, y, x) ;
    gphiu = np.ones((time,z_a,jpj,jpi),np.double);
    gphiu[0,0,:,:]=(NCin_dmn.variables['gphiu'][0,:,lon_cut:]).copy().astype(np.double);
    
#    double gphiv(time, z_a, y, x) ;
    gphiv = np.ones((time,z_a,jpj,jpi),np.double);
    gphiv[0,0,:,:]=(NCin_dmn.variables['gphiv'][0,:,lon_cut:]).copy().astype(np.double);
        
#    double gphif(time, z_a, y, x) ;
    gphif = np.ones((time,z_a,jpj,jpi),np.double);
    gphif[0,0,:,:]=(NCin_dmn.variables['gphif'][0,:,lon_cut:]).copy().astype(np.double);
    
#    double ff(time, z_a, y, x) ;
    ff = np.ones((time,z_a,jpj,jpi),np.double);
    ff[0,0,:,:]=(NCin_dmn.variables['ff_f'][0,:,lon_cut:]).copy().astype(np.double);

#    float nav_lat(y, x) ;
    nav_lat = np.ones((jpj,jpi),float);
    nav_lat[:,:] = (NCin_dmn.variables['nav_lat'][:,lon_cut:]).copy().astype(float);
    
#    float nav_lev(z) ;
    nav_lev = np.ones((jpk,),float);
    nav_lev[:] = (NCin_dmn.variables['nav_lev'][:jpk]).copy().astype(float);

#    float nav_lon(y, x) ;
    nav_lon = np.ones((jpj,jpi),float);
    nav_lon[:,:] = (NCin_dmn.variables['nav_lon'][:,lon_cut:]).copy().astype(float);
    
    
#    double e1f(time, z_a, y, x) ;
    e1f = np.ones((time,z_a,jpj,jpi),np.double);
    e1f[0,0,:,:] =  (NCin_dmn.variables['e1f'][0,:,lon_cut:]).copy().astype(np.double)

#    double e1t(time, z_a, y, x) ;
    e1t = np.ones((time,z_a,jpj,jpi),np.double);
    e1t[0,0,:,:] =  (NCin_dmn.variables['e1t'][0,:,lon_cut:]).copy().astype(np.double)
    
#    double e1u(time, z_a, y, x) ;
    e1u = np.ones((time,z_a,jpj,jpi),np.double);
    e1u[0,0,:,:] =  (NCin_dmn.variables['e1u'][0,:,lon_cut:]).copy().astype(np.double)

#    double e1v(time, z_a, y, x) ;
    e1v = np.ones((time,z_a,jpj,jpi),np.double);
    e1v[0,0,:,:] =  (NCin_dmn.variables['e1v'][0,:,lon_cut:]).copy().astype(np.double)
    
    
#    double e2f(time, z_a, y, x) ;
    e2f = np.ones((time,z_a,jpj,jpi),np.double);
    e2f[0,0,:,:] =  (NCin_dmn.variables['e2f'][0,:,lon_cut:]).copy().astype(np.double)
    
    
#    double e2t(time, z_a, y, x) ;
    e2t = np.ones((time,z_a,jpj,jpi),np.double)
    e2t[0,0,:,:] =  (NCin_dmn.variables['e2t'][0,:,lon_cut:]).copy().astype(np.double)
    
#    double e2u(time, z_a, y, x) ;
    e2u = np.ones((time,z_a,jpj,jpi),np.double);
    e2u[0,0,:,:] =  (NCin_dmn.variables['e2u'][0,:,lon_cut:]).copy().astype(np.double)
    
#    double e2v(time, z_a, y, x) ;
    e2v = np.ones((time,z_a,jpj,jpi),np.double);
    e2v[0,0,:,:] =  (NCin_dmn.variables['e2v'][0,:,lon_cut:]).copy().astype(np.double)
    

#    double e3t(time, z, y, x) ;
    e3t = np.ones((time,jpk,jpj,jpi),np.double)
    e3t[0,:,:,:] =  (NCin_dmn.variables['e3t_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)

#    double e3u(time, z, y, x) ;
    e3u = np.ones((time,jpk,jpj,jpi),np.double);
    e3u[0,:,:,:] =  (NCin_dmn.variables['e3u_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)

#    double e3v(time, z, y, x) ;
    e3v = np.ones((time,jpk,jpj,jpi),np.double);
    e3v[0,:,:,:] =  (NCin_dmn.variables['e3v_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)

#    double e3w(time, z, y, x) ;
    e3w = np.ones((time,jpk,jpj,jpi),np.double);
    e3w[0,:,:,:] =  (NCin_dmn.variables['e3w_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)


#    double tmask(time, z, y, x) ;
    tmask = np.ones((time,jpk,jpj,jpi),np.double);
    tmask[0,:,:,:]  = (NCin_msk.variables['tmask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    tmask[0,:,0, :] =0.;
    tmask[0,:,:, 0] =0.;
    tmask[0,:,jpj-1,:] =0.;
    tmask[0,:,:,jpi-1] =0.;
    tmask[0,jpk-1,:,:] =0.;

    if biscay_land:
        Lon=glamt[0,0,:,:];
        Lat=gphit[0,0,:,:];
        for j in range(jpj):
            for i in range(jpi):
                if(Lon[j,i]<0.0 and Lat[j,i]>42   ): tmask[0,:,j,i]=0.
                if(Lon[j,i]<-6. and Lat[j,i]>37.25): tmask[0,:,j,i]=0.

#    double umask(time, z, y, x) ;
    umask = np.ones((time,jpk,jpj,jpi),np.double);
    umask[0,:,:,:]  = (NCin_msk.variables['umask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    umask[0,jpk-1,:,:] =0.;

    if biscay_land:
        Lon=glamu[0,0,:,:];
        Lat=gphiu[0,0,:,:];
        for j in range(jpj):
            for i in range(jpi):
                if(Lon[j,i]<0.0 and Lat[j,i]>42   ): umask[0,:,j,i]=0.
                if(Lon[j,i]<-6. and Lat[j,i]>37.25): umask[0,:,j,i]=0.

#    double vmask(time, z, y, x) ;
    vmask = np.ones((time,jpk,jpj,jpi),np.double);
    vmask[0,:,:,:]  = (NCin_msk.variables['vmask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    vmask[0,jpk-1,:,:] =0.;

    if biscay_land:
        Lon=glamv[0,0,:,:];
        Lat=gphiv[0,0,:,:];
        for j in range(jpj):
            for i in range(jpi):
                if(Lon[j,i]<0.0 and Lat[j,i]>42   ): vmask[0,:,j,i]=0.
                if(Lon[j,i]<-6. and Lat[j,i]>37.25): vmask[0,:,j,i]=0.

    # mask the north sea
    maxlat = 47.28 # some latitude in the middle of France
    Lat = gphit[0,0,:,:];
    for j in range(jpj):
        for i in range(jpi):
            if (Lat[j,i]>maxlat): tmask[0,:,j,i]=0.

    Lat = gphiv[0,0,:,:];
    for j in range(jpj):
        for i in range(jpi):
            if (Lat[j,i]>maxlat): vmask[0,:,j,i]=0.

    Lat = gphiu[0,0,:,:];
    for j in range(jpj):
        for i in range(jpi):
            if (Lat[j,i]>maxlat): umask[0,:,j,i]=0.


    ##############################################################
    # write meshmask netcdf file !
    ##############################################################
    print(f"saving: {outfile}")
    ncOUT=NC.Dataset(outfile,"w");

    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',time)
    
    ncOUT.createDimension('x_a',x_a);
    ncOUT.createDimension('y_a',y_a);
    ncOUT.createDimension('z_a',z_a);

    ncvar    = ncOUT.createVariable('coastp','d',('y','x'))                  ; ncvar[:] = coastp;  
    ncvar    = ncOUT.createVariable('e1f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e1f   ;
    ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e1t   ;
    ncvar    = ncOUT.createVariable('e1u'   ,'d',('time','z_a', 'y', 'x')   ); ncvar[:] = e1u   ;
    ncvar    = ncOUT.createVariable('e1v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e1v   ;
    ncvar    = ncOUT.createVariable('e2f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e2f   ;
    ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e2t   ;
    ncvar    = ncOUT.createVariable('e2u'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2u   ;
    ncvar    = ncOUT.createVariable('e2v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2v   ;     

 
    ncvar    = ncOUT.createVariable('e3t_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3t   ;
    ncvar    = ncOUT.createVariable('e3u_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3u   ;
    ncvar    = ncOUT.createVariable('e3v_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3v   ;
    ncvar    = ncOUT.createVariable('e3w_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = e3w   ;   

    ncvar    = ncOUT.createVariable('ff'    ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = ff    ;      
    ncvar    = ncOUT.createVariable('fmask' ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = fmask ;    
    ncvar    = ncOUT.createVariable('gdept' ,'d',('time','z', 'y', 'x')) ; ncvar[:] = gdept ;
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

    ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = tmask 
    ncvar    = ncOUT.createVariable('umask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = umask 
    ncvar    = ncOUT.createVariable('vmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = vmask 
    ncOUT.close()
    

    
    
    
