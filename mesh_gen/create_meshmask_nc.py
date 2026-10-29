import numpy as np
import netCDF4 as NC
from flood_fill_ocean import flood_fill_nemo

# Script to create OGSTM meshmask file,
# from the original NEMO one.
# Specifically tailored (i.e. with hard-coded parameters)
# for med24v1 NEMO [curvilinear] regional grid, with
# a large portion of the North East Atlantic and 
# the Black Sea.
# Requirement: 
# NEMO meshmask must have full coordinate variables (glamt, gphit, etc.)
# i.e. no zeros or mask on land only mpi domains.

class OrigMask():

    def __init__(self, mask_file):
        NCin=NC.Dataset(mask_file,"r")
        self.Lon = (NCin.variables['glamt'][0,:,:]).copy()
        self.Lat = (NCin.variables['gphit'][0,:,:]).copy()
        self.nav_lev = (NCin.variables['nav_lev'][:]).copy().astype(float);
        self.nc_handler = NCin


    def get_wp(self):
        '''
        Prints some information about the number of waterpooints in longitude
        The aim is to provide some information about where to cut this mesh to get the ogstm meshmask.
        '''
        NCin=self.nc_handler
        jpi=NCin.dimensions['x'].size
        jpj=NCin.dimensions['y'].size
        tmask= (NCin.variables['tmask'][0,:,:,:]).copy().astype(np.int64) #shape is (141, 715, 1275)
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

    def getCutLocation(self, lon=None, lat=None):
        dist = 0
        if lon is not None:
            dist = dist + (self.Lon - lon)**2
        if lat is not None:
            dist = dist + (self.Lat - lat)**2
        j, i = np.unravel_index(np.nanargmin(dist), dist.shape)
        return j, i

def get_masks(NCin, Coords, jpk, jpj, jpi, Cut_Loc, i0, j0, blocked_edges):
    '''
    reads original tuvf-mask from NCin,
    Cooords are pairs of ((glam_, gphi_), ...), in order tuvw
    masks black sea with flood fill (args: i0, j0, blocked_edges),
    trims depth and lon,
    puts a 1-masked-cell frame all around the edges
    returns:
        tuvf-mask: np.array
    '''
    # read original masks
    orig_masks = [(NCin.variables[_mask][:]).copy().astype(np.double) for _mask in ('tmask', 'umask', 'vmask', 'fmask')]
    # flood-fill black sea
    masks = flood_fill_nemo(*orig_masks, i0, j0, blocked_edges=blocked_edges)
    # trim masks
    ilon_min, ilon_max, jlat_min, jlat_max, kdepth_max = Cut_Loc.values()
    masks = [_mask[:, :kdepth_max, jlat_min:jlat_max, ilon_min:ilon_max] for _mask in masks]
    #masks = [_mask[:,:jpk,:,lon_cut:] for _mask in masks]
    # frame masks
    for _mask in masks:
#        _mask[0,:,0, :] =0.;
#        _mask[0,:,:, 0] =0.;
#        _mask[0,:,jpj-1,:] =0.;
#        _mask[0,:,:,jpi-1] =0.;
#        _mask[0,jpk-1,:,:] =0.;
        _mask[0,:,0, :] = 0.
        _mask[0,:,:, 0] = 0.
        _mask[0,:,-1,:] = 0.
        _mask[0,:,:,-1] = 0.
        _mask[0,-1,:,:] = 0.
    # mask biscay and north sea (is not optional!)
    for (glam_, gphi_), _mask in zip(Coords, masks):
        Lon = glam_[0,0,:,:];
        Lat = gphi_[0,0,:,:];
        region = ((Lon < 0.0) & (Lat > 42)) | ((Lon < -6.0) & (Lat > 37.25))# | (Lat > 47.28)
        _mask[0, :, :, :] = np.where(region[None, :, :], 0.0, _mask[0, :, :, :])
    # output
    return masks

time = 1;
z_a  = 1;
y_a  = 1;
x_a  = 1;

def trim_var(vname, NCin, Cut_Loc, out_ndim):
    ilon_min, ilon_max, jlat_min, jlat_max, kdepth_max = Cut_Loc.values()
    X0 = (NCin.variables[vname][:]).copy().astype(np.double)
    print(f"{vname}.shape: {X0.shape}")
    if X0.ndim == 4: 
        # Dimensions: (time, depth, lat, lon)
        if X0.shape[1] == 1:
            # depth singleton
            kdepth_max = 1
        X = X0[0, :kdepth_max, jlat_min:jlat_max, ilon_min:ilon_max]
    elif X0.ndim == 3:
        # Dimensions: (time, lat, lon)
        X = X0[0, jlat_min:jlat_max, ilon_min:ilon_max]
    elif X0.ndim == 2:
        # Dimensions: (lat, lon)
        X = X0[jlat_min:jlat_max, ilon_min:ilon_max]
    elif X0.ndim == 1:
        # Dimensions: (depth)
        X = X0[:kdepth_max]
    else:
        raise ValueError(f"Unsupported dimensions for variable '{vname}': {X0.shape}")
    # add heading singleton dimensions
    X = X.reshape((1,) * (out_ndim - X.ndim) + X.shape)
    return X

def create_meshmask_nc(OrigMaskobj, outfile, Cut_Loc, depth_cut, biscay_land=True):
    
    NCin=OrigMaskobj.nc_handler

    ilon_min, ilon_max, jlat_min, jlat_max, kdepth_max = Cut_Loc.values()

    jpi = ilon_max - ilon_min
    jpj = jlat_max - jlat_min
    jpk = kdepth_max

#    jpi=NCin.dimensions['x'].size-lon_cut
#    jpj=NCin.dimensions['y'].size
#    jpk=NCin.dimensions['nav_lev'].size-depth_cut #116

    coastp = np.zeros((jpj,jpi),np.double)

##    double fmask(time, z, y, x) ;
#    fmask = np.ones((time,jpk,jpj,jpi),np.double);
#    fmask[0,:,:,:]=(NCin.variables['fmask'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
#
##    fmask[0,:,1, 0]    = 0.;
##    fmask[0,:,0, 1]    = 0.;
##    fmask[0,:,-1,0]    = 0.;
##    fmask[0,:,0,-1]    = 0.;

#    double gdept(time, z, y, x) ;
    #gdept              = np.zeros((time,jpk,jpj,jpi),np.double);
    #gdept[0,0:jpk,:,:] = (NCin.variables['gdept_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double);
    #gdept[0,:,:,:] = (NCin.variables['gdept_0'][0,:kdepth_max,jlat_min:jlat_max,ilon_min:ilon_max]).copy().astype(np.double);
    #gdept = (NCin.variables['gdept_0'][0,:kdepth_max,jlat_min:jlat_max,ilon_min:ilon_max]).copy().astype(np.double);
    gdept = trim_var('gdept_0', NCin, Cut_Loc, 4)

#    double gdepw(time, z, y_a, x_a) ;
    gdepw              = np.zeros((time,jpk,y_a,x_a),np.double);
    gdepw[0,0:jpk,0,0] = (NCin.variables['gdepw_1d'][0,:jpk]).copy().astype(np.double); 

#    double glamt(time, z_a, y, x) ;
    #glamt = np.ones((time,z_a,jpj,jpi),np.double);
    #glamt[0,0,:,:]=(NCin.variables['glamt'][0,:,lon_cut:]).copy().astype(np.double);
    glamt = trim_var('glamt', NCin, Cut_Loc, 4)
        
#    double glamu(time, z_a, y, x) ;
    #glamu = np.ones((time,z_a,jpj,jpi),np.double);
    #glamu[0,0,:,:]=(NCin.variables['glamu'][0,:,lon_cut:]).copy().astype(np.double);
    glamu = trim_var('glamu', NCin, Cut_Loc, 4)
        
#    double glamv(time, z_a, y, x) ;
    #glamv = np.ones((time,z_a,jpj,jpi),np.double);
    #glamv[0,0,:,:]=(NCin.variables['glamv'][0,:,lon_cut:]).copy().astype(np.double);
    glamv = trim_var('glamv', NCin, Cut_Loc, 4)

#    double glamf(time, z_a, y, x) ;
    #glamf = np.ones((time,z_a,jpj,jpi),np.double);
    #glamf[0,0,:,:]=(NCin.variables['glamf'][0,:,lon_cut:]).copy().astype(np.double);
    glamf = trim_var('glamf', NCin, Cut_Loc, 4)
        
#    double gphit(time, z_a, y, x) ;
    #gphit = np.ones((time,z_a,jpj,jpi),np.double);
    #gphit[0,0,:,:]=(NCin.variables['gphit'][0,:,lon_cut:]).copy().astype(np.double);
    gphit = trim_var('gphit', NCin, Cut_Loc, 4)
        
#    double gphiu(time, z_a, y, x) ;
    #gphiu = np.ones((time,z_a,jpj,jpi),np.double);
    #gphiu[0,0,:,:]=(NCin.variables['gphiu'][0,:,lon_cut:]).copy().astype(np.double);
    gphiu = trim_var('gphiu', NCin, Cut_Loc, 4)
    
#    double gphiv(time, z_a, y, x) ;
    #gphiv = np.ones((time,z_a,jpj,jpi),np.double);
    #gphiv[0,0,:,:]=(NCin.variables['gphiv'][0,:,lon_cut:]).copy().astype(np.double);
    gphiv = trim_var('gphiv', NCin, Cut_Loc, 4)
        
#    double gphif(time, z_a, y, x) ;
    #gphif = np.ones((time,z_a,jpj,jpi),np.double);
    #gphif[0,0,:,:]=(NCin.variables['gphif'][0,:,lon_cut:]).copy().astype(np.double);
    gphif = trim_var('gphif', NCin, Cut_Loc, 4)
    
#    double ff(time, z_a, y, x) ;
    #ff = np.ones((time,z_a,jpj,jpi),np.double);
    #ff[0,0,:,:]=(NCin.variables['ff_f'][0,:,lon_cut:]).copy().astype(np.double);
    ff = trim_var('ff_f', NCin, Cut_Loc, 4)

#    float nav_lat(y, x) ;
    #nav_lat = np.ones((jpj,jpi),float);
    #nav_lat[:,:] = (NCin.variables['nav_lat'][:,lon_cut:]).copy().astype(float);
    nav_lat = trim_var('nav_lat', NCin, Cut_Loc, 2)
    
#    float nav_lev(z) ;
    #nav_lev = np.ones((jpk,),float);
    #nav_lev[:] = (NCin.variables['nav_lev'][:jpk]).copy().astype(float);
    nav_lev = trim_var('nav_lev', NCin, Cut_Loc, 1)

#    float nav_lon(y, x) ;
    #nav_lon = np.ones((jpj,jpi),float);
    #nav_lon[:,:] = (NCin.variables['nav_lon'][:,lon_cut:]).copy().astype(float);
    nav_lon = trim_var('nav_lon', NCin, Cut_Loc, 2)
    
#    double e1f(time, z_a, y, x) ;
    #e1f = np.ones((time,z_a,jpj,jpi),np.double);
    #e1f[0,0,:,:] =  (NCin.variables['e1f'][0,:,lon_cut:]).copy().astype(np.double)
    e1f = trim_var('e1f', NCin, Cut_Loc, 4)

#    double e1t(time, z_a, y, x) ;
    #e1t = np.ones((time,z_a,jpj,jpi),np.double);
    #e1t[0,0,:,:] =  (NCin.variables['e1t'][0,:,lon_cut:]).copy().astype(np.double)
    e1t = trim_var('e1t', NCin, Cut_Loc, 4)
    
#    double e1u(time, z_a, y, x) ;
    #e1u = np.ones((time,z_a,jpj,jpi),np.double);
    #e1u[0,0,:,:] =  (NCin.variables['e1u'][0,:,lon_cut:]).copy().astype(np.double)
    e1u = trim_var('e1u', NCin, Cut_Loc, 4)

#    double e1v(time, z_a, y, x) ;
    #e1v = np.ones((time,z_a,jpj,jpi),np.double);
    #e1v[0,0,:,:] =  (NCin.variables['e1v'][0,:,lon_cut:]).copy().astype(np.double)
    e1v = trim_var('e1v', NCin, Cut_Loc, 4)
    
    
#    double e2f(time, z_a, y, x) ;
    #e2f = np.ones((time,z_a,jpj,jpi),np.double);
    #e2f[0,0,:,:] =  (NCin.variables['e2f'][0,:,lon_cut:]).copy().astype(np.double)
    e2f = trim_var('e2f', NCin, Cut_Loc, 4)
    
    
#    double e2t(time, z_a, y, x) ;
    #e2t = np.ones((time,z_a,jpj,jpi),np.double)
    #e2t[0,0,:,:] =  (NCin.variables['e2t'][0,:,lon_cut:]).copy().astype(np.double)
    e2t = trim_var('e2t', NCin, Cut_Loc, 4)
    
#    double e2u(time, z_a, y, x) ;
    #e2u = np.ones((time,z_a,jpj,jpi),np.double);
    #e2u[0,0,:,:] =  (NCin.variables['e2u'][0,:,lon_cut:]).copy().astype(np.double)
    e2u = trim_var('e2u', NCin, Cut_Loc, 4)
    
#    double e2v(time, z_a, y, x) ;
    #e2v = np.ones((time,z_a,jpj,jpi),np.double);
    #e2v[0,0,:,:] =  (NCin.variables['e2v'][0,:,lon_cut:]).copy().astype(np.double)
    e2v = trim_var('e2v', NCin, Cut_Loc, 4)
    

#    double e3t(time, z, y, x) ;
    #e3t = np.ones((time,jpk,jpj,jpi),np.double)
    #e3t[0,:,:,:] =  (NCin.variables['e3t_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    e3t = trim_var('e3t_0', NCin, Cut_Loc, 4)

#    double e3u(time, z, y, x) ;
    #e3u = np.ones((time,jpk,jpj,jpi),np.double);
    #e3u[0,:,:,:] =  (NCin.variables['e3u_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    e3u = trim_var('e3u_0', NCin, Cut_Loc, 4)

#    double e3v(time, z, y, x) ;
    #e3v = np.ones((time,jpk,jpj,jpi),np.double);
    #e3v[0,:,:,:] =  (NCin.variables['e3v_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    e3v= trim_var('e3v_0', NCin, Cut_Loc, 4)

#    double e3w(time, z, y, x) ;
    #e3w = np.ones((time,jpk,jpj,jpi),np.double);
    #e3w[0,:,:,:] =  (NCin.variables['e3w_0'][0,:jpk,:,lon_cut:]).copy().astype(np.double)
    e3w = trim_var('e3w_0', NCin, Cut_Loc, 4)


    # parameters to fill in all of marmara and black sea, 
    # and leave 1pt at dardanelles for boundary
    i0 = 980
    j0 = 320
    blocked_edges = [((299, 943), (299, 944))]

    # get masks
    #tmask = get__mask('tmask', NCin, glamt, gphit, jpk, jpj, jpi, lon_cut, i0, j0, blocked_edges)
    #umask = get__mask('umask', NCin, glamu, gphiu, jpk, jpj, jpi, lon_cut, i0, j0, blocked_edges)
    #vmask = get__mask('vmask', NCin, glamv, gphiv, jpk, jpj, jpi, lon_cut, i0, j0, blocked_edges)
    Coords = ((glamt, gphit), 
              (glamu, gphiu),
              (glamv, gphiv),
              (glamf, gphif))
    tmask, umask, vmask, fmask = get_masks(NCin, Coords, jpk, jpj, jpi, Cut_Loc, i0, j0, blocked_edges)

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
    

if __name__=="__main__":
    mesh_in = '/home/ggalli/Documents/RIVIERADE/PREPROC/med24v1_mesh_mask_BGC.nc'
    mesh_out = '/home/ggalli/Documents/RIVIERADE/PREPROC/med24v1_domain_cfg.nc'
    #mesh_in = '/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc'
    #mesh_out = '/leonardo_work/RIVRD_OGS/med24v1_domain_cfg.nc'

    M = OrigMask(mesh_in)
    jpi_orig = M.Lon.shape[1] 
    jpj_orig = M.Lat.shape[0]
    jpk_orig = M.nav_lev.shape[0]
    print(f"orig jpk, jpj, jpi: {jpk_orig} {jpj_orig} {jpi_orig}")

    # initialise mesh cut locations
    Cut_Loc = {'ilon_min': 0, 'ilon_max': jpi_orig,
               'jlat_min': 0, 'jlat_max': jpj_orig,
               'kdepth_max': jpk_orig}

    # cut out gibraltar
    lat_cut, lon_cut = M.getCutLocation(-8.7, 36.909)
    Cut_Loc['ilon_min'] = lon_cut
    # cut east of levantine
    lat_cut, lon_cut = M.getCutLocation(33.82, None)
    Cut_Loc['ilon_max'] = lon_cut
    # cut north of Med
    lat_cut, lon_cut = M.getCutLocation(None, 47.14)
    Cut_Loc['jlat_max'] = lat_cut
    # depth_cut #this has a different logic (-16)
    kdepth_cut = 16 
    Cut_Loc['kdepth_max'] = kdepth_cut

    Biscay_land = True

    #create_meshmask_nc(M, outfile, lon_cut, depth_cut, Biscay_land)
    create_meshmask_nc(M, mesh_out, Cut_Loc, Biscay_land)


    
    
