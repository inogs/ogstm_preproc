import numpy as np
import netCDF4 as nc
import logging
from commons.submask import SubMask
from basins import V2

class atmosphere():

    def __init__(self,mask,conf):
        '''
        Arguments:
        * mask * is a commons.Mask object
        * conf * is a configuration object, obtained by read_configure
        '''

        logging.info("Atmosphere start calculation")
        _,jpj,jpi = mask.shape

        self.eas = self.eas + self.aeg + self.adn + self.ads;
        sup_mask = mask.cut_at_level(0)
        tmask = mask.mask_at_level(0)
        wes = SubMask(V2.eas,maskobject=sup_mask).mask_at_level(0)
        eas = SubMask(V2.eas,maskobject=sup_mask).mask_at_level(0)
        
        aux01=wes.copy()
        aux02=eas.copy()
        for jj in range(1,jpj-1):
            for ji in range (1,jpi-1):
                if ((wes[jj,ji] == 0) or (eas[jj,ji] == 0) ) and (tmask[jj,ji] == 1):
                    if ((wes[jj+1,ji] ) or (wes[jj-1,ji] ) or (wes[jj,ji+1]) or (wes[jj,ji-1]) ):
                        aux01[jj,ji] = 1;
                    else:
                        aux02[0,jj,ji] = 1;


        wes = aux01;
        eas = aux02;

        Nwes = 0;
        Neas = 0;
        Cell_Area = mask.area
        e3t       = mask.dz
        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                Nwes = Nwes + Cell_Area[jj,ji]*e3t[jj,ji]*wes[jj,ji]
                Neas = Neas + Cell_Area[jj,ji]*e3t[jj,ji]*eas[jj,ji]
                #Nwes = Nwes + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,jj,ji]*self.wes[0,jj,ji];
                #Neas = Neas + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,jj,ji]*self.eas[0,jj,ji];


        self.nitrate   = np.zeros((jpj,jpi),np.float32);
        self.phosphate = np.zeros((jpj,jpi),np.float32);

        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                if (wes[jj,ji]):
                    self.nitrate[jj,ji] = conf.n3n_wes/Nwes
                    self.phosphate[jj,ji] = conf.po4_wes/Nwes
                if (eas[jj,ji]):
                    self.nitrate[jj,ji] = conf.n3n_eas/Neas
                    self.phosphate[jj,ji] = conf.po4_eas/Neas

        logging.info("Atmosphere finish calculation")

    def write_atm_netcdf(self,mask,outdir):
        _,jpj,jpi = mask.shape

        l_tmask = ~ mask.mask_at_level(0)

        area = mask.area()
        w = 1.0e+09;
        t = 1/(365 * 86400);
        n = 1/14
        p = 1/31
        totP = 0;
        totN = 0;
        totN_KTy =0;
        totP_KTy =0;
        e3t = mask.dz
        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                VolCell1=area[jj,ji]*e3t[jj,ji]#self._mesh_father.e3t[0,0,jj,ji];
                totN = totN + self.nitrate[jj,ji]*VolCell1;
                totP = totP + self.phosphate[jj,ji]*VolCell1;
                totN_KTy = totN_KTy + self.nitrate[jj,ji]*(1.e-3/n)*VolCell1;
                totP_KTy = totP_KTy + self.phosphate[jj,ji]*(1.e-3/p)*VolCell1;
                cn = w*t;
                cp = w*t;
                self.nitrate[jj,ji]   = self.nitrate[jj,ji]*cn;
                self.phosphate[jj,ji] = self.phosphate[jj,ji]*cp;


        fileOUT = outdir + "/ATM_yyyy0630-00:00:00.nc"
        self.nitrate[l_tmask] = np.nan
        self.phosphate[l_tmask] = np.nan

        ncfile = nc.Dataset(fileOUT, "w", format="NETCDF4")
        ncfile.createDimension('lon', jpi)
        ncfile.createDimension('lat', jpj)
        n = ncfile.createVariable('atm_N3n','f', ('lat','lon'))
        n[:] = self.nitrate[:]
        p = ncfile.createVariable('atm_N1p','f', ('lat','lon'))
        p[:] = self.phosphate[:]
        setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
        setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
        setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
        setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
        ncfile.close()

        logging.info("Atmosphere Netcdf writed")
