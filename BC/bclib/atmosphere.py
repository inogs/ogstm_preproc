import numpy as np
import netCDF4 as nc
import logging
from bitsea.commons.submask import SubMask
from bitsea.basins import V2
from bitsea.basins import V2 as OGS
from bitsea.basins.basin import ComposedBasin


class atmosphere:

    def __init__(self, mask, conf):
        """
        Arguments:
        * mask * is a `bitsea.commons.Mask` object
        * conf * is a configuration object, obtained by importing `config.py`

        nitrate and phosphate are fields expressed in Mmol/m3/y
        """

        logging.info("Atmosphere start calculation")
        _, jpj, jpi = mask.shape

        EAS = ComposedBasin('eas', [V2.eas, V2.adr1, V2.adr2, V2.aeg], 'East with marginal seas')

        # Here we use basin object features to run faster
        EAS.region.border_longitudes = [9, 36]
        EAS.region.border_latitudes = [30, 46]
        V2.wes.region.border_longitudes = [-6, 18]
        V2.wes.region.border_latitudes = [32, 45]
        V2.med.region.border_longitudes = [-6, 36]
        V2.med.region.border_latitudes = [30, 46]
        # ---------------------------------------------
        print("wes")
        wes = SubMask(V2.wes, mask)
        print("eas")
        eas = SubMask(EAS, mask)
        print("done")

        e3t = mask.dz
        Volume = mask.area * e3t[0]
        Nwes = Volume[wes[0]].sum()
        Neas = Volume[eas[0]].sum()

        self.nitrate = np.zeros((jpj, jpi), np.float32)
        self.phosphate = np.zeros((jpj, jpi), np.float32)

        self.nitrate[wes[0]] = conf.n3n_wes / Nwes
        self.phosphate[wes[0]] = conf.po4_wes / Nwes
        self.nitrate[eas[0]] = conf.n3n_eas / Neas
        self.phosphate[eas[0]] = conf.po4_eas / Neas

        logging.info("Atmosphere finish calculation")

    def check(self, mask):

        _, jpj, jpi = mask.shape
        n = 1. / 14
        p = 1. / 31
        vol_layer0 = mask.area * mask.dz[0]
        totP = np.sum(self.phosphate * vol_layer0)
        totN = np.sum(self.nitrate * vol_layer0)
        totN_KTy = totN * (1.e-3 / n)
        totP_KTy = totP * (1.e-3 / p)

        return totN, totP, totN_KTy, totP_KTy

    def write_netcdf(self, mask, outdir, concentration_dimension="Volume"):
        '''
        Arguments
        * mask                     * is the Mask object
        * outdir                   * (string) where ATM_yyyy0630-00:00:00.nc will be written.
        *  concentration_dimension * string can be "Volume" or "Area"
        '''

        assert concentration_dimension == "Volume" or concentration_dimension == "Area"
        # Check is done using volume concentrations
        totN, totP, totN_KTy, totP_KTy = self.check(mask)
        # conversion from Mmol/y to mmol/s
        w = 1.0e+09
        t = 1. / (365 * 86400)
        cn = w * t
        self.nitrate = self.nitrate * cn
        self.phosphate = self.phosphate * cn

        if concentration_dimension == "Volume":
            units = "mmol/m3"
            hunits = "nmol/m3"
            comment = "transport model can use these atmospherical contributions as they are"
        elif concentration_dimension == "Area":
            units = "mmol/m2"
            comment = "transport model can use thes atmospherical contributions by dividing them by the depth of the first cell"
            self.nitrate = self.nitrate * mask.dz[0]
            self.phosphate = self.phosphate * mask.dz[0]
        else:
            raise ValueError(
                "concentration_dimension must be 'Volume' or 'Area'"
            )

        # now is ready to dump

        fileOUT = outdir + "ATM_yyyy0630-00:00:00.nc"
        _, jpj, jpi = mask.shape

        ncfile = nc.Dataset(fileOUT, "w", format="NETCDF4")
        ncfile.createDimension('lon', jpi)
        ncfile.createDimension('lat', jpj)
        n = ncfile.createVariable('atm_N3n', 'f', ('lat', 'lon'))
        n[:] = self.nitrate[:]
        setattr(n, "units", units)
        p = ncfile.createVariable('atm_N1p', 'f', ('lat', 'lon'))
        p[:] = self.phosphate[:]
        setattr(p, "units", units)
        setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
        setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
        setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
        setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
        setattr(ncfile, 'comments', comment)
        ncfile.close()

        logging.info("Atmosphere Netcdf written.")

class atmosphere_seas:

    def __init__(self, mask, conf):
        """
        Arguments:
        * mask * is a `bitsea.commons.Mask` object
        * conf * is a configuration object, obtained by importing `config.py`

        nitrate and phosphate are fields read in mmol m-2 s-1
        """

        logging.info("Atmosphere start calculation")        
        jpk, jpj, jpi = mask.shape
        Mask0 = mask.cut_at_level(0)
        mask200= mask.mask_at_level(200)

        mydtype=[('sub','U20')]
        VARS=['NOw','NOs','NCw','NCs','POw','POs','PCw','PCs']
        mydtype.extend( [(s,np.float32) for s in VARS] )
        A=np.loadtxt(conf.file_atm,dtype=mydtype, skiprows=2)

        A[A['sub']=='alb']['NCw']

        dtype = [(sub.name, bool) for sub in OGS.Pred]
        SUB = np.zeros((jpj,jpi),dtype=dtype)
        for sub in OGS.Pred:
            sbmask         = SubMask(sub, mask)
            SUB[sub.name]  = sbmask[0,:,:]

        self.N3n_win = np.zeros((jpj,jpi),np.float32)
        self.N3n_sum = np.zeros((jpj,jpi),np.float32)
        self.N1p_win = np.zeros((jpj,jpi),np.float32)
        self.N1p_sum = np.zeros((jpj,jpi),np.float32)


        for sub in OGS.Pred:
            if sub.name=='atl':continue
            ii_sub = A['sub']==sub.name
            mask = SUB[sub.name]
            self.N3n_win[mask] = A[ii_sub]['NCw']
            self.N3n_sum[mask] = A[ii_sub]['NCs']
            self.N1p_win[mask] = A[ii_sub]['PCw']
            self.N1p_sum[mask] = A[ii_sub]['PCs']

            mask = SUB[sub.name] & mask200
            self.N3n_win[mask] = A[ii_sub]['NOw']
            self.N3n_sum[mask] = A[ii_sub]['NOs']
            self.N1p_win[mask] = A[ii_sub]['POw']
            self.N1p_sum[mask] = A[ii_sub]['POs']


        logging.info("Atmosphere finish calculation")

    def check(self, mask, season="summer"):

        _, jpj, jpi = mask.shape
        n = 1. / 14
        p = 1. / 31
        vol_layer0 = mask.area * mask.dz[0]
        if season=="summer":
            totP = np.sum(self.N1p_sum * vol_layer0)
            totN = np.sum(self.N3n_sum * vol_layer0)
        else:
            totP = np.sum(self.N1p_win * vol_layer0)
            totN = np.sum(self.N3n_win * vol_layer0)

        totN_KTy = totN * (1.e-3 / n)
        totP_KTy = totP * (1.e-3 / p)

        return totN, totP, totN_KTy, totP_KTy
    


    def write_netcdf(self, mask, outdir):
        '''
        Arguments
        * mask                     * is the Mask object
        * outdir                   * (string) where ATM_yyyy0630-00:00:00.nc will be written.
        *  concentration_dimension * string can be "Volume" or "Area"
        '''


        def dump_netcdf(fileOUT, N3n, N1p, totN_KTy, totP_KTy, totP, totN, totNcomment):
            print(fileOUT,flush=True)
            with nc.Dataset(fileOUT, "w", format="NETCDF4") as ncfile:
                ncfile.createDimension('lon', jpi)
                ncfile.createDimension('lat', jpj)
                n = ncfile.createVariable('atm_N3n', 'f', ('lat', 'lon'))
                n[:] = N3n[:]
                setattr(n, "units", "mmol/m2/s")
                p = ncfile.createVariable('atm_N1p', 'f', ('lat', 'lon'))
                p[:] = N1p[:]
                setattr(p, "units", "mmol/m2/s")
                setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
                setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
                setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
                setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
                setattr(ncfile, 'comments', totNcomment)

            
        _, jpj, jpi = mask.shape        
        fileOUT = outdir + "ATM_yyyy0331-00:00:00.nc"
        dump_netcdf(fileOUT, self.N3n_win, self.N1p_win, 0,0,0,0, "winter values")
        fileOUT = outdir + "ATM_yyyy0401-00:00:00.nc"
        dump_netcdf(fileOUT, self.N3n_sum, self.N1p_sum, 0,0,0,0, "summer values")
        fileOUT = outdir + "ATM_yyyy0930-00:00:00.nc"
        dump_netcdf(fileOUT, self.N3n_sum, self.N1p_sum, 0,0,0,0, "summer values")
        fileOUT = outdir + "ATM_yyyy1001-00:00:00.nc"
        dump_netcdf(fileOUT, self.N3n_win, self.N1p_win, 0,0,0,0, "winter values")
 

        logging.info("Atmosphere Netcdf written.")