import numpy as np
import numpy.matlib as npmat
import netCDF4 as nc
import logging
import code
import matplotlib.pyplot as plt


class sub_mesh:
    """
        Class sub mesh

    """

    def __init__(self,mesh,ncfile):
        self.path = ncfile
        self._mesh_father = mesh
        self._extract_information()
        logging.info("submesh builded")

    def _extract_information(self):
        try:
            self.ncfile = nc.Dataset(self.path, 'r')
        except:
            print("SUNMASK NOT FOUND")
            exit()
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()

    def atmosphere(self):


        logging.info("Atmosphere start calculation")
        jpk = self._mesh_father.tmask_dimension[1]
        jpj = self._mesh_father.tmask_dimension[2]
        jpi = self._mesh_father.tmask_dimension[3]

        self.eas = self.eas + self.aeg + self.adn + self.ads;

        aux01=self.wes;
        aux02=self.eas;
        for jj in range(1,jpj-2):
            for ji in range (1,jpi-2):
                if ((self.wes[0,jj,ji] == 0) or (self.eas[0,jj,ji] == 0) ) and (self._mesh_father.tmask[0,0,jj,ji] == 1):
                    if ((self.wes[0,jj+1,ji] == 1) or (self.wes[0,jj-1,ji] == 1) or (self.wes[0,jj,ji+1] == 1) or (self.wes[0,jj,ji-1] == 1) ):
                        self.wes[0,jj,ji] = 1;
                    else:
                        aux02[0,jj,ji] = 1;


        self.wes = aux01;
        self.eas = aux02;

        Nwes = 0;
        Neas = 0;
        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                Nwes = Nwes + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,jj,ji]*self.wes[0,jj,ji];
                Neas = Neas + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,jj,ji]*self.eas[0,jj,ji];

        lon = self._mesh_father.nav_lon
        lat = self._mesh_father.nav_lat
        self.atm = np.zeros((jpj,jpi,2));
        a = self._mesh_father.input_data


        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                if (self.wes[0,jj,ji] == 1):
                    self.atm[jj,ji,0] = a.n3n_wes/Nwes
                    self.atm[jj,ji,1] = a.po4_wes/Nwes
                if (self.eas[0,jj,ji] == 1):
                    self.atm[jj,ji,0] = a.n3n_eas/Neas
                    self.atm[jj,ji,1] = a.po4_eas/Neas

        self.write_atm_netcdf()

        logging.info("Atmosphere finish calculation")

    def write_atm_netcdf(self):

        jpk = self._mesh_father.tmask_dimension[1]
        jpj = self._mesh_father.tmask_dimension[2]
        jpi = self._mesh_father.tmask_dimension[3]
        l_tmask = ~ self._mesh_father.tmask[0,0][:].astype(np.bool)

        ntra_atm_a  = self.atm[:,:,0];
        phos_atm_a  = self.atm[:,:,1];
        area = self._mesh_father.e1t[0,0,:,:]*self._mesh_father.e2t[0,0,:,:]
        w = 1.0e+09;
        t = 1/(365 * 86400);
        n = 1/14
        p = 1/31
        totP = 0;
        totN = 0;
        totN_KTy =0;
        totP_KTy =0;

        for jj in range(0,jpj-1):
             for ji in range(0,jpi-1):
                    VolCell1=area[jj,ji]*self._mesh_father.e3t[0,0,jj,ji];
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
