import numpy as np
import numpy.matlib as npmat
import netCDF4 as nc
import logging
import code
import matplotlib.pyplot as plt

class bmesh:

    """
        Class bounmesh

    """

    def __init__(self,mesh,ncfile):
        self.path = ncfile
        self._mesh_father = mesh
        logging.info("bounmesh builded")

    def load_bounmask(self):

        try:
            bncfile = nc.Dataset(self.path, 'r')
            print("LOADING BMASK FROM FILE")
            for i in bncfile.dimensions:
                print(bncfile.dimensions[i])
                setattr(self, bncfile.dimensions[i].name, bncfile.dimensions[i].size)
            for i in bncfile.variables:
                b = bncfile.variables[i][:].copy()
                setattr(self, i, b)
                print(i)
            bncfile.close()
        except:
            print("BOUNMASK NOT FOUND")
            exit()

    def write_netcdf(self):
        logging.info("Start bounmesh nc file write")
        time = 1
        ncfile = nc.Dataset(self.path, 'w')
        ncfile.createDimension('x',self._mesh_father.x)
        ncfile.createDimension('y',self._mesh_father.y)
        ncfile.createDimension('z',self._mesh_father.z)
        ncfile.createDimension('time',time)
        ncfile.createDimension('waterpoints',self.water_points)
        ncfile.createDimension('dim3',3)

        navlon_wnc = ncfile.createVariable('nav_lon', 'f', ('y','x'))
        navlon_wnc[:] = self._mesh_father.nav_lon[:]
        navlat_wnc = ncfile.createVariable('nav_lat', 'f', ('y','x'))
        navlat_wnc[:] = self._mesh_father.nav_lat[:]
        navlev_wnc = ncfile.createVariable('nav_lev', 'f', 'z')
        navlev_wnc[:] = self._mesh_father.nav_lev[:]
        for i in range(0,self.nudg):
            corrf =[1.,1.,1.,1.,1.01,1.01,1.];
            aux= self.resto[i,:,:,:]*corrf[i];
            np.transpose(aux, (2, 1, 0)).shape
            name = "re" + self.vnudg[i][0]
            #print(name)
            resto_wnc = ncfile.createVariable(name, 'f4', ('time','z','y','x'))
            resto_wnc[0,:] = aux
        idx_inv_wnc = ncfile.createVariable('index', 'i', ('time','z','y','x'))
        idx_inv_wnc[0,:] = self.idx # self.idx[:]
        idx_rev_wnc = ncfile.createVariable('index_inv', 'i', ('waterpoints', 'dim3'))
        idx_l_inv = self.idx_inv +1
        idx_rev_wnc[:] = idx_l_inv
        ncfile.close()
        logging.info("bounmesh nc file writed")
