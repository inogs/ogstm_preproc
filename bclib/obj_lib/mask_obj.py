
import numpy as np
from scipy.io import netcdf as nc
from bclib.io_lib import excel_obj as xlsobj
import logging

class lateral_bc:
    
    def __init__(self,file_nutrients):
        self.path = file_nutrients
        self._extract_information()
        logging.info("lateral_bc builded") 
    
    def _extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()
    
    
        
    
class river_data:
    
    def __init__(self,mesh,file_river,file_runoff):
        self.path_river = file_river
        self.path_runoff = file_runoff
        self._mesh_father = mesh 
        self._extract_information()
        logging.info("river_data builded") 
        
    def _extract_information(self):
        logging.debug(self.path_river)
        logging.debug(self.path_runoff)
        work_sheet  = self._mesh_father.input_data.river_data_sheet
        
        range_montly = range(7,19)
        range_coord  = [1,2]
        
        river_excel_file = xlsobj.xlsx(self.path_river)
        self.river_montly_mod = river_excel_file.read_spreadsheet_allrow("monthly",range_montly)
        self.river_coordr = river_excel_file.read_spreadsheet_allrow("monthly",range_coord)
        self.nrivers = len(self.river_coordr[:])
        self.river_collected_data = {}
        for data_t in self._mesh_father.input_data.river_data_sheet:
            river_sheet_collected_data = {}
            x_range_coord  = [1]
            y_range = range(7,48)
            self.river_years = river_excel_file.read_spreadsheet_range(data_t,x_range_coord,y_range,"i")
            count = 0
            x_range = range(2,39)
            ry = river_excel_file.read_spreadsheet_range(data_t,x_range,y_range)
            for y in self.river_years[0][:]:
                river_sheet_collected_data[str(y)] = ry[:,count].copy()
                logging.debug(str(y))
            self.river_collected_data[data_t] =  river_sheet_collected_data.copy()    
        logging.debug("End river data collection")
              
        roundoff_excel_file = xlsobj.xlsx(self.path_runoff)
        self.roundoff_montly_mod = river_excel_file.read_spreadsheet_allrow("monthly",range_montly)
        self.roundoff_coordr = river_excel_file.read_spreadsheet_allrow("monthly",range_coord)
        self.nRoundoff = len(self.river_coordr[:])
        self.roundoff_collected_data = {}
        for data_t in self._mesh_father.input_data.river_data_sheet:
            roundoff_sheet_collected_data = {}
            x_range_coord  = [1]
            y_range = range(7,48)
            self.roundoff_years = river_excel_file.read_spreadsheet_range(data_t,x_range_coord,y_range,"i")
            count = 0
            x_range = range(2,39)
            ry = roundoff_excel_file.read_spreadsheet_range(data_t,x_range,y_range)
            for y in self.river_years[0][:]:
                roundoff_sheet_collected_data[str(y)] = ry[:,count].copy()
            self.roundoff_collected_data[data_t] = roundoff_sheet_collected_data.copy()
        logging.debug("End roundoff data collection")
            
    def _coast_line_mask(self, mask):
        [rows, cols] = np.nonzero(mask)
        # [rows, cols] = find(mask);
        nSea = len(rows);
        coast=bool(mask);
        
        for i in range(1,nSea):
            row=rows[i]
            col=cols[i]        
            irows=[-1, 0, 1] + row
            icols=[-1, 0, 1] + col
            localmask = mask[irows,icols]
            
            coast[row,col] = sum(localmask[:]) < 9
        return coast    
                   
    def map_contribute_on_sea(self):
        mask1 = self._mesh_father.tmask[0,1,:,:]
        mask2 = self._mesh_father.tmask[0,2,:,:]
        coast = self._coast_line_mask(mask1) and self._coast_line_mask(mask2)
        
        loncm = self._mesh_father.nav_lon(coast)
        latcm = self._mesh_father.nav_lat(coast)
        
        [coastline_row_index, coastline_col_index] = np.nonzero(coast)
        georef4 = np.array(coastline_row_index, coastline_col_index, loncm, latcm) 
        
        data_types = self._mesh_father.input_data.river_data_sheet
        n_data_types = len(data_types)
        
        georef = np.zeros(self.nRivers,5)
        for jr in range (1,self.nRivers):
            lon_river = self.river_coordr(jr,1)
            lat_river = self.river_coordr(jr,2)
            dist = (loncm-lon_river)**2 + (latcm-lat_river)**2
            ind = np.amin(dist)
            # w = np.min(dist)
            georef[jr,:]=np.array(jr,georef4[ind,:])
        river_georef = georef
        
        m=np.zeros(self.nRivers,12)
        self.river_data={}
        
        for data_type in data_types :
            years_data={}
            for ic in self.river_years :
                for r in range (1,self.nRivers) :
                    ry = self.river.river_collected_data["data_type"]["ic"][r]
                    m[r,:] =  (self.river_montly_mod[r,:]/100)*12*ry
                years_data[str(ic)]=m.copy()
            river_data[data_type]=years_data.copy() 
        
        n_coast_cells = len(loncm)
        
        indexes = np.zeros(n_coast_cells,1)
        georef = np.zeros(self.nRivers,5)
        
        m=np.zeros(self.nRivers,12)
        self.river_data={}
        
        for data_type in data_types :
            years_data={}
            for ic in self.river_years :
                for r in range (1,self.nRivers) :
                    ry = self.river.river_collected_data["data_type"]["ic"][r]
                    m[r,:] =  (self.river_montly_mod[r,:]/100)*12*ry
                years_data[str(ic)]=m.copy()
            river_data[data_type]=years_data.copy()                    
        
        
    










    


class boun_mesh:

    """
        Class bounmesh

    """

    def __init__(self,mesh,ncfile):
        self.path = ncfile
        self._mesh_father = mesh 
        logging.info("bounmesh builded") 

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
        logging.info("bounmesh nc file writed") 
        
     
    


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
        logging.info("atmosphere finish calculation") 
        

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
        self.gibilterra = lateral_bc(self.input_data.file_nutrients)
        self.river = river_data(self,self.input_data.file_river,self.input_data.file_runoff)
        logging.info("mesh builded") 
        

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
        logging.info("bounmesh generation ended") 
        
        
        
    def bc(self):
        pass
        
        
        
        
        
        
        
        
        
