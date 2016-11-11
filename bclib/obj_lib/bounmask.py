import numpy as np
import netCDF4 as nc
import logging

class bounmask():

    def __init__(self,conf):
        self.config  = conf
        self.resto   = None
        self.idx     = None
        self.idx_inv = None
        logging.info("bounmask builded")
        

    def generate(self,mask):

        """ Generation of bounmask fields
        resto,idx, and idx_inv
         """


        vnudg = self.config.variables
        rdpmin = np.float64(self.config.rdpmin)
        rdpmax = np.float64(self.config.rdpmax)
        glamt  = np.float64(mask.xlevels)
        nudg = len(vnudg)
        jpk, jpjglo, jpiglo = mask.shape


        resto = np.zeros((nudg,jpk,jpjglo,jpiglo), dtype=np.float64);
        for jn in range(nudg):
            xlim = vnudg[jn][1]
            ii = glamt<xlim
            resto[jn,:,ii] = 1./(rdpmin*86400.)
        for jn in range(nudg):
            xlim = vnudg[jn][1]
            ii = (glamt > xlim)  & (glamt <= self.config.end_nudging)
            reltim = rdpmin + (rdpmax-rdpmin)*(glamt-xlim)/(self.config.end_nudging-xlim);
            for jk in range(jpk):
                resto[jn,jk,ii] = 1./(reltim[ii]*86400.)

        resto[:,~mask.mask] = 1.e+20;
        count = 0
        idx = np.zeros((jpk,jpjglo,jpiglo),dtype=np.int)
        self.water_points = np.sum(mask.mask)
        idx_inv = np.zeros((self.water_points,3),dtype=np.int);

        for jk in range(jpk):
            for jj in range(jpjglo):
                for ji in range(jpiglo):
                    if mask.mask[jk,jj,ji]:
                        idx[jk,jj,ji] = count+1;
                        idx_inv[count,0]=jk+1;
                        idx_inv[count,1]=jj+1;
                        idx_inv[count,2]=ji+1;
                        count=count+1;


        idx[~mask.mask] = 0;
        self.idx = idx
        self.idx_inv = idx_inv
        self.resto = resto

        #self.write_netcdf()
        logging.info("bounmesh generation ended")



    def load(self,varname):
        '''

        Arguments:
        * varname * (string) of a variable in bounmask.nc
        If varname is a "resto" variable or "index"
        the returned value is a 3d array
        If varname is "inde_inv" the returned value is a 2d array

        Returns:
        numpy array
        '''

        bncfile = nc.Dataset(self.config.file_bmask, 'r')
        print("LOADING from bounmask")
        if (varname =="index_inv"):
            M = np.array(bncfile[varname])
        else:
            M = np.array(bncfile[varname][0,:,:,])
        bncfile.close()
        return M


    def write_netcdf(self,mask):
        logging.info("Start bounmesh nc file write")
        vnudg = self.config.variables
        nudg = len(vnudg)
        ncfile = nc.Dataset(self.config.file_bmask, 'w')
        jpk,jpj,jpi = mask.shape
        ncfile.createDimension('x',jpi)
        ncfile.createDimension('y',jpj)
        ncfile.createDimension('z',jpk)
        ncfile.createDimension('time',1)
        ncfile.createDimension('waterpoints',self.water_points)
        ncfile.createDimension('dim3',3)

        navlon_wnc = ncfile.createVariable('nav_lon', 'f', ('y','x'))
        navlon_wnc[:] = mask.xlevels
        navlat_wnc = ncfile.createVariable('nav_lat', 'f', ('y','x'))
        navlat_wnc[:] = mask.ylevels
        navlev_wnc = ncfile.createVariable('nav_lev', 'f', 'z')
        navlev_wnc[:] = mask.zlevels
        for jn in range(nudg):
            corrf =[1.,1.,1.,1.,1.01,1.01,1.];
            aux= self.resto[jn,:,:,:]*corrf[jn];
            np.transpose(aux, (2, 1, 0)).shape
            name = "re" + vnudg[jn][0]
            resto_wnc = ncfile.createVariable(name, 'f4', ('time','z','y','x'))
            resto_wnc[0,:] = aux
            setattr(resto_wnc,'missing_value',1.e+20)
        idx_inv_wnc = ncfile.createVariable('index', 'i', ('time','z','y','x'))
        idx_inv_wnc[0,:] = self.idx
        setattr(idx_inv_wnc,'missing_value',0)
        idx_rev_wnc = ncfile.createVariable('index_inv', 'i', ('waterpoints', 'dim3'))
        idx_l_inv = self.idx_inv +1
        idx_rev_wnc[:] = idx_l_inv
        ncfile.close()
        logging.info("bounmask.nc file writed")

if __name__ == '__main__':
    from bclib.io_lib  import read_configure
    from commons.mask import Mask
    conf = read_configure.elaboration(json_input="../../conf24.json")
    conf.file_mask="../../meshmask_872.nc"
    conf.file_bmask="bounmask.nc"
    conf.active_bmask = True
    TheMask = Mask(conf.file_mask)
    B=bounmask(conf)
    B.generate(TheMask)
    B.write_netcdf(TheMask)
    
    B_old = bounmask(conf)
    B_old.config.file_bmask = "/Users/gbolzon/Documents/workspace/ogs_bounday_conditions/bounmask8.nc"
    index= B_old.load('index')
    restoN1p=B_old.load('reN1p')
    d=B.resto[0,:,:,:]-restoN1p

