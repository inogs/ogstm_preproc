import numpy as np
import netCDF4 as nc
import logging
from commons.submask import SubMask
from basins import V2
from basins.basin import ComposedBasin

class atmosphere():

    def __init__(self,mask,conf):
        '''
        Arguments:
        * mask * is a commons.Mask object
        * conf * is a configuration object, obtained by importing config.py

        nitrate and phoshate are fields expressed in Mmol/m3/y
        '''

        logging.info("Atmosphere start calculation")
        _,jpj,jpi = mask.shape

        mask.cut_at_level(0)
        EAS=ComposedBasin('eas',[V2.eas,V2.adr1,V2.adr2,V2.aeg],'East with marginal seas')
        wes = SubMask(V2.wes,maskobject=mask).mask_at_level(0)
        eas = SubMask(EAS,   maskobject=mask).mask_at_level(0)
        
        e3t       = mask.dz
        Volume = mask.area*e3t[0]
        Nwes = Volume[wes].sum()
        Neas = Volume[eas].sum()

        self.nitrate   = np.zeros((jpj,jpi),np.float32);
        self.phosphate = np.zeros((jpj,jpi),np.float32);

        self.nitrate[wes]  = conf.n3n_wes/Nwes
        self.phosphate[wes]= conf.po4_wes/Nwes
        self.nitrate[eas]  = conf.n3n_eas/Neas
        self.phosphate[eas]= conf.po4_eas/Neas

#         for jj in range(0,jpj-1):
#             for ji in range(0,jpi-1):
#                 if (wes[jj,ji]):
#                     self.nitrate[jj,ji] = conf.n3n_wes/Nwes
#                     self.phosphate[jj,ji] = conf.po4_wes/Nwes
#                 if (eas[jj,ji]):
#                     self.nitrate[jj,ji] = conf.n3n_eas/Neas
#                     self.phosphate[jj,ji] = conf.po4_eas/Neas

        logging.info("Atmosphere finish calculation")

    def check(self,mask):

        _,jpj,jpi = mask.shape
        n = 1./14
        p = 1./31
        totP     = .0
        totN     = .0
        totN_KTy = .0
        totP_KTy = .0
        for jj in range(jpj):
            for ji in range(jpi):
                VolCell1=mask.area[jj,ji]*mask.dz[0]#[jj,ji]#self._mesh_father.e3t[0,0,jj,ji];
                totN = totN + self.nitrate[jj,ji]  *VolCell1;
                totP = totP + self.phosphate[jj,ji]*VolCell1;
                totN_KTy = totN_KTy + self.nitrate[jj,ji]  *(1.e-3/n)*VolCell1;
                totP_KTy = totP_KTy + self.phosphate[jj,ji]*(1.e-3/p)*VolCell1;
        return totN,totP, totN_KTy, totP_KTy



    def write_netcdf(self,mask,outdir):

        totN,totP, totN_KTy, totP_KTy = self.check(mask)
        #conversion from Mmol/y to mmol/s
        w = 1.0e+09
        t = 1./(365 * 86400)
        cn = w*t
        self.nitrate = self.nitrate*cn
        self.phosphate = self.phosphate*cn
        # now is ready to dump

        fileOUT = outdir + "/ATM_yyyy0630-00:00:00.nc"
        _,jpj,jpi = mask.shape

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
