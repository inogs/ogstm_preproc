import numpy as np
import numpy.matlib as npmat
import netCDF4 as nc
import logging
import code
import matplotlib.pyplot as plt

class bounmask():

    def __init__(self,ncfile, conf):
        self.path = ncfile
        logging.info("bounmesh builded")
        

    def generate_bounmask(self):

        """ This fuction generate bounmask """

        if self.input_data.active_bmask == True :
            #input var
            bm = self.bounmesh
            bm.vnudg = self.input_data.variables
            rdpmin = self.input_data.rdpmin
            rdpmax = self.input_data.rdpmax
            #print(type(self.y),type(self.x))
            self.glamt = self.glamt.reshape(int(self.y),int(self.x))
            bm.nudg = len(bm.vnudg)
            bm.jpk = self.tmask_dimension[1]
            bm.jpjglo = self.tmask_dimension[2]
            bm.jpiglo = self.tmask_dimension[3]
            #print(bm.nudg,bm.jpk,bm.jpjglo,bm.jpiglo)
            bm.resto = np.zeros((bm.nudg,bm.jpk,bm.jpjglo,bm.jpiglo));

            for jk in range(0,bm.jpk):

                for jn in range(0,bm.nudg):
                    for jj in range(0,bm.jpjglo):
                        for ji in range(0,bm.jpiglo):

                            if (self.glamt[jj][ji] < bm.vnudg[jn][1]):

                                bm.resto[jn,jk,jj,ji]=1./(rdpmin*86400.);

                for jn in range(0,bm.nudg):
                    for jj in range(0,bm.jpjglo):
                        for ji in range(0,bm.jpiglo):
                            if (self.glamt[jj][ji] > bm.vnudg[jn][1]) and (self.glamt[jj][ji] <= self.input_data.end_nudging):
                                reltim = rdpmin + (rdpmax-rdpmin)*(self.glamt[jj][ji]-bm.vnudg[jn][1])/(self.input_data.end_nudging-bm.vnudg[jn][1]);
                                bm.resto[jn,jk,jj,ji] = 1./(reltim*86400.);


            bm.resto[:,self.tmask[0] == 0] = 1.e+20;
            count = 0
            bm.idx = np.zeros((bm.jpk,bm.jpjglo,bm.jpiglo),dtype=np.int)
            bm.water_points = np.sum(self.tmask)
            bm.idx_inv = np.zeros((bm.water_points,3),dtype=np.int);

            for jk in range(bm.jpk):
                for jj in range(bm.jpjglo):
                    for ji in range(bm.jpiglo):
                        if self.tmask[0,jk,jj,ji] == 1.0:
                            bm.idx[jk,jj,ji] = count+1;
                            bm.idx_inv[count,0]=jk+1;
                            bm.idx_inv[count,1]=jj+1;
                            bm.idx_inv[count,2]=ji+1;
                            count=count+1;


            bm.idx[self.tmask[0] == 0] = 0;
            self.bounmesh.write_netcdf()
            logging.info("bounmesh generation ended")




        else :
            logging.info("bounmesh generation disabled")


    def load_bounmask(self):
        '''
        Returns only the index 3d integer matrix
        '''
        try:
            bncfile = nc.Dataset(self.path, 'r')
            print("LOADING BMASK FROM FILE")
            index = np.array(bncfile['index'][0,:,:,])
            bncfile.close()
            return index
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
