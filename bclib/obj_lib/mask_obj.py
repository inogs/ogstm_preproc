
import numpy as np
from scipy.io import netcdf as nc


class boun_mesh:

    """
        Class bounmesh

    """

    def __init__(self, ncfile):
        self.path = ncfile


    def write_netcdf(self):
        print("PATH = ", self.path)
        ncfile = nc.netcdf_file(self.path, 'w')
        ncfile.createDimension('x',self.x)
        ncfile.createDimension('y',self.y)
        ncfile.createDimension('z',self.z)
        ncfile.createDimension('time',self.time)
        ncfile.createDimension('waterpoints',self.water_points)
        ncfile.createDimension('dim3',3)
        ncfile.close()
# dimensions:
#     x = 362 ;
#     y = 128 ;
#     z = 72 ;
#     time = 1 ;
#     waterpoints = 782035 ;
#     dim3 = 3 ;

# variables:
#     float nav_lon(y, x) ;
#     float nav_lat(y, x) ;
#     float nav_lev(z) ;
#     double reN1p(time, z, y, x) ;
#     double reN3n(time, z, y, x) ;
#     double reO2o(time, z, y, x) ;
#     double reN5s(time, z, y, x) ;
#     double reO3c(time, z, y, x) ;
#     double reO3h(time, z, y, x) ;
#     double reN6r(time, z, y, x) ;
#     int index(time, z, y, x) ;
#     int index_inv(waterpoints, dim3) ;


class sub_mesh:
    """
        Class sub mesh

    """

    def __init__(self, ncfile):
        self.path = ncfile
        self.extract_information()


    def extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            setattr(self, i, self.ncfile.variables[i][:])


class mesh:
    """
        Class mesh

    """

    def __init__(self, ncf_mesh,ncf_submesh,ncf_bounmesh):
        self.path = ncf_mesh
        self.extract_information()
        self.submesh = sub_mesh(ncf_submesh)
        self.bounmesh = boun_mesh(ncf_bounmesh)

    def extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            setattr(self, i, self.ncfile.variables[i][:])

    def generate_boundmask(self,elab):

        """

        This fuction generate boundmask
        from elab we use :
        vnudg  : array of variable
        end_nudging : elab.end_nudging
        rdmpmin : elab.rdmpmin
        rdmpmax : elab.rdmpmin

        """
        
        # aggiungere parametri di default per rdmpmin = 1/24.
        # e rdmpmax = 90.;
        
        #input var
        bm = self.bounmesh
        bm.x = self.x
        bm.y = self.y
        bm.z = self.z
        bm.time = self.time
        
        vnudg = elab.variables
        rdpmin = elab.rdpmin
        rdpmax = elab.rdpmax
        self.glamt = self.glamt.reshape(self.y,self.x)
        nudg = len(vnudg)
        tmask_dimension = self.tmask.shape
        bm.jpk = tmask_dimension[1]
        bm.jpjglo = tmask_dimension[2]
        bm.jpiglo = tmask_dimension[3]
        
        bm.resto = np.zeros((nudg,bm.jpk,bm.jpjglo,bm.jpiglo));
        
        for jk in range(1,bm.jpk):
            
            for jn in range(0,nudg):
                for jj in range(0,bm.jpjglo):
                    for ji in range(0,bm.jpiglo): 
                                       
                        if (self.glamt[jj][ji] < vnudg[jn][1]):
                            
                            bm.resto[jn,jk,jj,ji]=1./(rdpmin*86400.);
        
            for jn in range(0,nudg):
                for jj in range(0,bm.jpjglo):
                    for ji in range(0,bm.jpiglo):
                        if (self.glamt[jj][ji] > vnudg[jn][1]) and (self.glamt[jj][ji] <= elab.end_nudging):
                            reltim = rdpmin + (rdpmax-rdpmin)*(self.glamt[jj][ji]-vnudg[jn][1])/(elab.end_nudging-vnudg[jn][1]);
                            bm.resto[jn,jk,jj,ji] = 1./(reltim*86400.);
        
        
        bm.resto[:,self.tmask[0] == 0] = 1.e+20;
        count = 0
        bm.idx = np.zeros((bm.jpk,bm.jpjglo,bm.jpiglo),dtype=np.int)
        bm.water_points = np.sum(self.tmask)
        bm.idx_inv = np.zeros((bm.water_points,3),dtype=np.int);
        
        for jk in range(1,bm.jpk):
            for jj in range(1,bm.jpjglo):
                for ji in range(1,bm.jpiglo):
                    if self.tmask[0,jk,jj,ji] == 1:
                        count=count+1;
                        bm.idx[jk,jj,ji] = count;
                        bm.idx_inv[count,:]=[jk,jj,ji];
        

        bm.idx[self.tmask[0] == 0] = 0;
        
        
        
        
        
        
        
        
        
        
