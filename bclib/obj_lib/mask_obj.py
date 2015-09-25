
import numpy as np
from scipy.io import netcdf as nc
from bclib.io_lib import excel_obj as xlsobj


class gibilterra:
    def __init__(self):
        pass
    
class river_data:
    
    def __init__(self,mesh,file_river,file_runoff):
        self.path_river = file_river
        self.path_runoff = file_runoff
        self._mesh_father = mesh 
        self._extract_information()
        
    def _extract_information(self):
        print(self.path_river)
        print(self.path_runoff)
        work_sheet  = self._mesh_father.input_data.river_data_sheet
        
        range_montly = range(7,19)
        range_coord  = [1,2]
        
        river_excel_file = xlsobj.xlsx(self.path_river)
        montly_mod = river_excel_file.read_spreadsheet_allrow("monthly",range_montly)
        coordr = river_excel_file.read_spreadsheet_allrow("monthly",range_coord)
        range_coord  = [1]
        a = river_excel_file.read_spreadsheet_allcols("monthly",range_coord)
        #roundoff_excel_file = xlsobj.xlsx(self.path_river)
        
        
        
    
    
    


class boun_mesh:

    """
        Class bounmesh

    """

    def __init__(self,mesh,ncfile):
        self.path = ncfile
        self._mesh_father = mesh 

    def write_netcdf(self):
        
        ncfile = nc.netcdf_file(self.path, 'w')
        ncfile.createDimension('x',self._mesh_father.x)
        ncfile.createDimension('y',self._mesh_father.y)
        ncfile.createDimension('z',self._mesh_father.z)
        ncfile.createDimension('time',self._mesh_father.time)
        ncfile.createDimension('waterpoints',self.water_points)
        ncfile.createDimension('dim3',3)
        
        navlon_wnc = ncfile.createVariable('nav_lon', 'f', ('y','x'))
        navlon_wnc = self._mesh_father.nav_lon
        navlat_wnc = ncfile.createVariable('nav_lat', 'f', ('y','x'))
        navlat_wnc = self._mesh_father.nav_lat
        navlev_wnc = ncfile.createVariable('nav_lev', 'f', 'z')
        navlev_wnc = self._mesh_father.nav_lev
        for i in range(0,self.nudg):
            corrf =[1.,1.,1.,1.,1.01,1.01,1.];
            aux= self.resto[i,:,:,:]*corrf[i];
            np.transpose(aux, (2, 1, 0)).shape
            name = "re" + self.vnudg[i][0]
            resto_wnc = ncfile.createVariable(name, 'd', ('time','z','y','x'))
            resto_wnc = aux
        idx_inv_wnc = ncfile.createVariable('index', 'i', ('time','z','y','x'))
        idx_inv_wnc = self.idx_inv
        idx_rev_wnc = ncfile.createVariable('index_inv', 'i', ('waterpoints', 'dim3'))
        idx_rev_wnc = np.transpose(self.idx, (2, 1, 0)).shape  
        ncfile.close()
        
     
    


class sub_mesh:
    """
        Class sub mesh

    """

    def __init__(self,mesh,ncfile):
        self.path = ncfile
        self._mesh_father = mesh
        self._extract_information()


    def _extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()
        
    def atmosphere(self):
        
        
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
                        self.wes[1,jj,ji] = 1;
                    else:
                        aux02[1,jj,ji] = 1;
        
        self.wes = aux01;
        self.eas = aux02;
        
        Nwes = 0;
        Neas = 0;
        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                Nwes = Nwes + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,0,0]*self.wes[0,jj,ji];
                Neas = Neas + self._mesh_father.e1t[0,0,jj,ji]*self._mesh_father.e2t[0,0,jj,ji]*self._mesh_father.e3t[0,0,0,0]*self.eas[0,jj,ji];
        
        count = 0; 
        idx = np.zeros((jpk,jpj,jpi));
        for jk in range(0,jpk-1):
            for jj in range(0,jpj-1):
                for ji in range(0,jpi-1):
                    if (self._mesh_father.tmask[0,jk,jj,ji] == 1):
                        count=count+1;
                        idx[jk,jj,ji] = count;
        
        counter =0 ; 
        for jj in range(1,jpj):
            for ji in range(1,jpi):
                if (self.wes[1,jj,ji] == 1):
                    counter=counter+1;
                if (self.eas[1,jj,ji] == 1):
                    counter=counter+1;
        
        
        self.atm = np.zeros((counter,7));
        lon = self._mesh_father.nav_lon
        lat = self._mesh_father.nav_lat
        a = self._mesh_father.input_data
        counter = 0
        for jj in range(0,jpj-1):
            for ji in range(0,jpi-1):
                jr=idx[0,jj,ji];
                if (self.wes[0,jj,ji] == 1):
                    counter=counter+1;
                    self.atm[counter,:] = [jr,ji,jj,lon[jj,ji],lat[jj,ji],a.n3n_wes/Nwes,a.po4_wes/Nwes]; 
                if (self.eas[0,jj,ji] == 1):
                    counter=counter+1 ;
                    self.atm[counter,:]=[jr,ji,jj,lon[jj,ji],lat[jj,ji],a.n3n_eas/Neas,a.po4_eas/Neas]; 

        

class mesh:
    """
        Class mesh
        
    """

    def __init__(self,input):# ncf_mesh,ncf_submesh,ncf_bounmesh):
        self.input_data = input
        self.path = self.input_data.file_mask
        self._extract_information()
        self.submesh = sub_mesh(self,self.input_data.file_submask)
        self.bounmesh = boun_mesh(self,self.input_data.file_bmask)
        self.river = river_data(self,self.input_data.file_river,self.input_data.file_runoff)

    def _extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()
        
        self.tmask_dimension = self.tmask.shape
        
        
    def generate_boundmask(self):

        """ This fuction generate boundmask """
        
        #input var
        bm = self.bounmesh        
        bm.vnudg = self.input_data.variables
        rdpmin = self.input_data.rdpmin
        rdpmax = self.input_data.rdpmax
        self.glamt = self.glamt.reshape(self.y,self.x)
        bm.nudg = len(bm.vnudg)
        bm.jpk = self.tmask_dimension[1]
        bm.jpjglo = self.tmask_dimension[2]
        bm.jpiglo = self.tmask_dimension[3]
        
        bm.resto = np.zeros((bm.nudg,bm.jpk,bm.jpjglo,bm.jpiglo));
        
        for jk in range(1,bm.jpk):
            
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
        
        for jk in range(1,bm.jpk):
            for jj in range(1,bm.jpjglo):
                for ji in range(1,bm.jpiglo):
                    if self.tmask[0,jk,jj,ji] == 1:
                        count=count+1;
                        bm.idx[jk,jj,ji] = count;
                        bm.idx_inv[count,:]=[jk,jj,ji];
        

        bm.idx[self.tmask[0] == 0] = 0;
        
        
        
    def bc(self):
        pass
        
        
        
        
        
        
        
        
        
