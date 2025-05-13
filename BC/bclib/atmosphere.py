import numpy as np
import netCDF4 as nc
import logging
from bitsea.commons.submask import SubMask
from bitsea.basins import V2
from bitsea.basins.basin import ComposedBasin

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
        V2.med.region.border_longitudes = [-6,36]
        V2.med.region.border_latitudes  = [30,46]
        #---------------------------------------------
        print("wes")
        wes = SubMask(V2.wes, mask0).mask_at_level(0)
        print("eas")
        eas = SubMask(EAS, mask0).mask_at_level(0)
        print("med")
        med = SubMask(V2.med,mask0).mask_at_level(0)
        print("done")
        
        e3t       = mask0.dz
        Volume = mask0.area*e3t[0]
        Nwes = Volume[wes].sum()
        Neas = Volume[eas].sum()
        Nmed = Volume[med].sum()

        self.nitrate   = np.zeros((jpj,jpi),np.float32);
        self.phosphate = np.zeros((jpj,jpi),np.float32);
        self.HgII      = np.zeros((jpj,jpi),np.float32);
        self.MMHg      = np.zeros((jpj,jpi),np.float32);
        self.Hg0atm    = np.zeros((jpj,jpi),np.float32);

        self.nitrate[wes]  = conf.n3n_wes/Nwes
        self.phosphate[wes]= conf.po4_wes/Nwes
        self.nitrate[eas]  = conf.n3n_eas/Neas
        self.phosphate[eas]= conf.po4_eas/Neas
        self.HgII[med]     = conf.HgII_med/Nmed
        self.MMHg[med]     = conf.MMHg_med/Nmed
        self.Hg0atm[med]  = conf.Hg0atm_med/Nmed
#       self.Hg0atm[eas]   = conf.Hg0atm_eas/Neas
#       self.HgII[eas]     = conf.HgII_eas/Neas
#       self.MMHg[eas]     = conf.MMHg_eas/Neas




        logging.info("Atmosphere finish calculation")

    def check(self,mask):

        _,jpj,jpi = mask.shape
        n = 1./14
        p = 1./31
        h=1./200.59
        totP     = .0
        totN     = .0
        totHgII  = .0
        totMMHg  = .0
        totN_KTy = .0
        totP_KTy = .0
        totHgII_KTy=.0
        totMMHg_KTy=.0

        for jj in range(jpj):
            for ji in range(jpi):
                VolCell1=mask.area[jj,ji]*mask.dz[0]
                totN = totN + self.nitrate[jj,ji]  *VolCell1
                totP = totP + self.phosphate[jj,ji]*VolCell1
                totHgII = totHgII + self.HgII[jj,ji]  *VolCell1
                totMMHg = totMMHg + self.MMHg[jj,ji]  *VolCell1

                totN_KTy = totN_KTy + self.nitrate[jj,ji]  *(1.e-3/n)*VolCell1
                totP_KTy = totP_KTy + self.phosphate[jj,ji]*(1.e-3/p)*VolCell1
                totHgII_KTy = totHgII_KTy + self.HgII[jj,ji]*(1.e-3/h)*VolCell1
                totMMHg_KTy = totMMHg_KTy + self.MMHg[jj,ji]*(1.e-3/h)*VolCell1
        return totN,totP,totHgII,totMMHg, totN_KTy, totP_KTy,totHgII_KTy, totMMHg_KTy



    def write_netcdf(self,mask,outdir,concentration_dimension="Volume"):
        '''
        Arguments
        * mask                     * is the Mask object
        * outdir                   * (string) where ATM_yyyy0630-00:00:00.nc will be written.
        *  concentration_dimension * string can be "Volume" or "Area"
        '''

        assert concentration_dimension=="Volume" or concentration_dimension =="Area"
        # Check is done using volume concentrations
        totN,totP,totHgII, totMMHg, totN_KTy, totP_KTy,totHgII_KTy, totMMHg_KTy = self.check(mask)
        #conversion from Mmol/y to mmol/s
        w = 1.0e+09
        w_hg = 1.0e+15
        t = 1./(365 * 86400)
        cn = w*t
        cn_hg=w_hg*t     #CHECK MERCURY CONVERSIONS
        self.nitrate   = self.nitrate*cn
        self.phosphate = self.phosphate*cn
        self.HgII      = self.HgII*cn_hg
        self.MMHg      = self.MMHg*cn_hg

        if concentration_dimension == "Volume":
            units   ="mmol/m3"
            hunits  ="nmol/m3"
            comment="transport model can use these atmospherical contributions as they are"
        if concentration_dimension =="Area":
            units = "mmol/m2"
            hunits = "nmol/m2"
            comment = "transport model can use thes atmospherical contributions by dividing them by the depth of the first cell"
            self.nitrate   = self.nitrate   * mask.dz[0]
            self.phosphate = self.phosphate * mask.dz[0]
            self.HgII = self.HgII * mask.dz[0]
            self.MMHg = self.MMHg * mask.dz[0]


        # now is ready to dump

        fileOUT = outdir + "ATM_yyyy0630-00:00:00_LOWHg.nc"
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
        h1 = ncfile.createVariable('atm_Hg2','f', ('lat','lon'))
        h1[:] = self.HgII[:]
        setattr(h1,"units", hunits)
        h2 = ncfile.createVariable('atm_MMHg','f', ('lat','lon'))
        h2[:] = self.MMHg[:]
        setattr(h2,"units", hunits)
        setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
        setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
        setattr(ncfile, 'ATM_HgII_MassBalance_kTON_y', totHgII_KTy)
        setattr(ncfile, 'ATM_MMHg_MassBalance_kTON_y', totMMHg_KTy)
        setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
        setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
        setattr(ncfile, 'ATM_HgII_MassBalance_Mmol_y', totHgII)
        setattr(ncfile, 'ATM_MMHg_MassBalance_Mmol_y', totMMHg)
        setattr(ncfile, 'comments', comment)
        ncfile.close()

        logging.info("Atmosphere Netcdf written.")

        fileOUT = outdir + "atm_Hg0_yyyy0107-00:00:00.nc"
        _,jpj,jpi = mask.shape

        ncfile = nc.Dataset(fileOUT, "w", format="NETCDF4")
        ncfile.createDimension('lon', jpi)
        ncfile.createDimension('lat', jpj)
        h = ncfile.createVariable('atm_Hg0','f', ('lat','lon'))
        h[:] = self.Hg0atm[:]
        setattr(h,"units", "ug m-3")
        ncfile.close()

        logging.info("Hg0 at mNetcdf written.")
