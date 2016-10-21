import numpy as np

import netCDF4 as nc
import logging
from commons.submask import SubMask
from basins import V2



class sub_mesh:
    """
        Class sub mesh

    """

    def __init__(self,mesh,ncfile):
        
        logging.info("submesh builded")

    def atmosphere(self,mask):


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


        self.atm = np.zeros((jpj,jpi,2));
        a = self._mesh_father.input_data


        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                if (wes[jj,ji]):
                    self.atm[jj,ji,0] = a.n3n_wes/Nwes
                    self.atm[jj,ji,1] = a.po4_wes/Nwes
                if (eas[jj,ji]):
                    self.atm[jj,ji,0] = a.n3n_eas/Neas
                    self.atm[jj,ji,1] = a.po4_eas/Neas

        self.write_atm_netcdf()

        logging.info("Atmosphere finish calculation")

    def write_atm_netcdf(self,mask):
        _,jpj,jpi = mask.shape

        l_tmask = ~ mask.mask_at_level(0)

        ntra_atm_a  = self.atm[:,:,0];
        phos_atm_a  = self.atm[:,:,1];
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
                    totN = totN + ntra_atm_a[jj,ji]*VolCell1;
                    totP = totP + phos_atm_a[jj,ji]*VolCell1;
                    totN_KTy = totN_KTy + ntra_atm_a[jj,ji]*(1.e-3/n)*VolCell1;
                    totP_KTy = totP_KTy + phos_atm_a[jj,ji]*(1.e-3/p)*VolCell1;
                    cn = w*t;
                    cp = w*t;
                    ntra_atm_a[jj,ji] = ntra_atm_a[jj,ji]*cn;
                    phos_atm_a[jj,ji] = phos_atm_a[jj,ji]*cp;

        for yCO2 in (range(self._mesh_father.input_data.simulation_start_time,
                            self._mesh_father.input_data.simulation_end_time)):

            fileOUT = self._mesh_father.input_data.dir_out + "/ATM_" + str(yCO2) + "0630-00:00:00.nc"
            #print(fileOUT)
            #map_co2 = np.dot(np.ones([self.input_data.jpj, self.input_data.jpi]), rcp85[count])
            ntra_atm_a[l_tmask] = np.nan
            phos_atm_a[l_tmask] = np.nan

            #ncfile = nc.netcdf_file(fileOUT, 'w')
            ncfile = nc.Dataset(fileOUT, "w", format="NETCDF4")
            ncfile.createDimension('lon', jpi)
            ncfile.createDimension('lat', jpj)
            n = ncfile.createVariable('atm_N3n','f', ('lat','lon'))
            n[:] = ntra_atm_a[:]
            p = ncfile.createVariable('atm_N1p','f', ('lat','lon'))
            p[:] = phos_atm_a[:]
            setattr(ncfile, 'ATM_P_MassBalance_kTON_y', totP_KTy)
            setattr(ncfile, 'ATM_N_MassBalance_kTON_y', totN_KTy)
            setattr(ncfile, 'ATM_P_MassBalance_Mmol_y', totP)
            setattr(ncfile, 'ATM_N_MassBalance_Mmol_y', totN)
            ncfile.close()

        logging.info("Atmosphere Netcdf writed")
