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

        mask0 = mask.cut_at_level(0)
        EAS=ComposedBasin('eas',[V2.eas,V2.adr1,V2.adr2,V2.aeg],'East with marginal seas')
        # Here we use basin objet features to run faster
        EAS.region.border_longitudes    = [ 9,36]
        EAS.region.border_latitudes     = [30,46]
        V2.wes.region.border_longitudes = [-6,18]
        V2.wes.region.border_latitudes  = [32,45]
        #---------------------------------------------
        print("wes")
        wes = SubMask(V2.wes,maskobject=mask0).mask_at_level(0)
        print("eas")
        eas = SubMask(EAS,   maskobject=mask0).mask_at_level(0)
        print("done")
        
        e3t       = mask0.dz
        Volume = mask0.area*e3t[0]
        Nwes = Volume[wes].sum()
        Neas = Volume[eas].sum()

        self.nitrate   = np.zeros((jpj,jpi),np.float32);
        self.phosphate = np.zeros((jpj,jpi),np.float32);

        self.nitrate[wes]  = conf.n3n_wes/Nwes
        self.phosphate[wes]= conf.po4_wes/Nwes
        self.nitrate[eas]  = conf.n3n_eas/Neas
        self.phosphate[eas]= conf.po4_eas/Neas


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
                VolCell1=mask.area[jj,ji]*mask.dz[0]
                totN = totN + self.nitrate[jj,ji]  *VolCell1
                totP = totP + self.phosphate[jj,ji]*VolCell1
                totN_KTy = totN_KTy + self.nitrate[jj,ji]  *(1.e-3/n)*VolCell1
                totP_KTy = totP_KTy + self.phosphate[jj,ji]*(1.e-3/p)*VolCell1
        return totN,totP, totN_KTy, totP_KTy



    def write_netcdf(self,mask,outdir,concentration_dimension="Volume"):
        '''
        Arguments
        * mask                     * is the Mask object
        * outdir                   * (string) where ATM_yyyy0630-00:00:00.nc will be written.
        *  concentration_dimension * string can be "Volume" or "Area"
        '''

        assert concentration_dimension=="Volume" or concentration_dimension =="Area"
        # Check is done using volume concentrations
        totN,totP, totN_KTy, totP_KTy = self.check(mask)
        #conversion from Mmol/y to mmol/s
        w = 1.0e+09
        t = 1./(365 * 86400)
        cn = w*t
        self.nitrate   = self.nitrate*cn
        self.phosphate = self.phosphate*cn

        if concentration_dimension == "Volume":
            units   ="mmol/m3"
            comment="transport model can use these atmospherical contributions as they are"
        if concentration_dimension =="Area":
            units = "mmol/m2"
            comment = "transport model can use thes atmospherical contributions by dividing them by the depth of the first cell"
            self.nitrate   = self.nitrate   * mask.dz[0]
            self.phosphate = self.phosphate * mask.dz[0]


        # now is ready to dump

        fileOUT = outdir + "ATM_yyyy0630-00:00:00.nc"
        _,jpj,jpi = mask.shape

        ncfile = nc.Dataset(fileOUT, "w", format="NETCDF4")
        ncfile.createDimension('lon', jpi)
        ncfile.createDimension('lat', jpj)
        n = ncfile.createVariable('atm_N3n','f', ('lat','lon'))
        n[:] = self.nitrate[:]
        setattr(n,"units", units)
        p = ncfile.createVariable('atm_N1p','f', ('lat','lon'))
        p[:] = self.phosphate[:]
        setattr(p,"units", units)
        setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
        setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
        setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
        setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
        setattr(ncfile, 'comments', comment)
        ncfile.close()

        logging.info("Atmosphere Netcdf written.")
