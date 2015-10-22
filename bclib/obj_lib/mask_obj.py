
import numpy as np
import numpy.matlib as npmat
from scipy.io import netcdf as nc
from scipy import interpolate as npint
from bclib.io_lib import excel_obj as xlsobj
import logging
from numpy import dtype

class lateral_bc:
    
    def __init__(self,mesh,file_nutrients):
        self.path = file_nutrients
        self._mesh_father = mesh 
        self._extract_information()
        self._convert_information()
        self.season = (["0215-12:00:00","0515-12:00:00",
                        "0815-12:00:00","1115-12:00:00"])
        logging.info("lateral_bc builded") 
    
    def _extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()
        
    def _convert_information(self):
        jpt = 4
        size_nutrients = np.zeros(4,dtype = np.int)
        size_nutrients[0] = jpt
        size_nutrients[1:] = self._mesh_father.tmask_dimension[1:]
        self.phos = np.zeros(size_nutrients)
        self.ntra = np.zeros(size_nutrients)
        self.dox = np.zeros(size_nutrients)
        self.sica = np.zeros(size_nutrients)
        self.dic = np.zeros(size_nutrients)
        self.alk = np.zeros(size_nutrients)
        tmask4d = np.zeros(size_nutrients)
        
        n_lev = self._mesh_father.tmask_dimension[1]
        vp_phos = np.zeros((4,n_lev)); 
        vp_ntra = np.zeros((4,n_lev)); 
        vp_sica = np.zeros((4,n_lev)); 
        vp_dox  = np.zeros((4,n_lev)); 
        vp_dic  = np.zeros((4,n_lev));
        vp_alk  = np.zeros((4,n_lev));
        
        nav_lev_in=np.concatenate(([0],self.lev1,[5500]))
        
        vp_phos_in=np.column_stack((self.N1p[:,1],self.N1p,self.N1p[:,len(self.lev1)-2]))
        
        vp_ntra_in=np.column_stack((self.N3n[:,1],self.N3n,self.N3n[:,len(self.lev1)-2]))
                             
        vp_sica_in=np.column_stack((self.N5s[:,1],self.N5s,self.N5s[:,len(self.lev1)-2]))
        
        vp_dox_in=np.column_stack((self.O2o[ :,1],self.O2o,self.O2o[:,len(self.lev1)-2]))
        end = len(self.ALK)
        
        nav_lev_in2 = np.append([0],self.lev2.T)
        nav_lev_in2 = np.append(nav_lev_in2,[5500])
        
        vp_dic_in = np.append(self.DIC[1],self.DIC[:].T)
        vp_dic_in = np.append(vp_dic_in,self.DIC[end-1]) 
        
        vp_alk_in = np.append(self.ALK[1],self.ALK[:].T)
        vp_alk_in = np.append(vp_alk_in,self.ALK[end-1]);
                       
        for i in range(4):
            jj= (~ np.isnan(vp_phos_in[i,:])) | ( vp_phos_in[i,:] < 1e19 ) ;
            vp_phos[i,:] = np.interp(self._mesh_father.nav_lev, nav_lev_in[jj], vp_phos_in[i,jj])
            jj=~ np.isnan(vp_ntra_in[i,:]) | ( vp_ntra_in[i,:] < 1e19 );  
            vp_ntra[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_ntra_in[i,jj]);
            jj=~ np.isnan(vp_dox_in[i,:]) | ( vp_dox_in[i,:] < 1e19 );
            vp_dox[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_dox_in[ i,jj]);
            jj=~ np.isnan(vp_sica_in[i,:]) | ( vp_sica_in[i,:] < 1e19);
            vp_sica[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_sica_in[i,jj]);
            jj=~ np.isnan(vp_dic_in[:]) | ( vp_dic_in[:] < 1e19);
            vp_dic[i,:]  = np.interp(self._mesh_father.nav_lev,nav_lev_in2[jj],vp_dic_in[jj]);
            jj=~ np.isnan(vp_alk_in[:]) | ( vp_alk_in[:] < 1e19) ; 
            vp_alk[i,:]  = np.interp(self._mesh_father.nav_lev,nav_lev_in2[jj],vp_alk_in[jj]);

# 
# %Loop on time seasonal
            for jt in range(jpt):
                for jk in range(n_lev):
                    self.phos[jt,jk][:,:] = vp_phos[jt,jk];
                    self.ntra[jt,jk,:,:] = vp_ntra[jt,jk];
                    self.dox[jt,jk,:,:] = vp_dox[jt,jk];
                    self.sica[jt,jk,:,:] = vp_sica[jt,jk];
                    self.dic[jt,jk,:,:]  = vp_dic[jt,jk];
                    self.alk[jt,jk,:,:]  = vp_alk[jt,jk];

 
            for jt in range(jpt):
                    tmask4d[jt,:,:,:]=self._mesh_father.tmask[:,:,:];
                    
                    
            self.phos[tmask4d == 0]=1.e+20;
            self.ntra[tmask4d == 0]=1.e+20;
            self.dox[tmask4d == 0]=1.e+20;
            self.sica[tmask4d == 0]=1.e+20;
            self.dic[tmask4d == 0]=1.e+20;
            self.alk[ tmask4d == 0]=1.e+20;

    
           
    
        
    
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
            self.river_years = self.river_years[0][:]
            count = 0
            x_range = range(2,39)
            ry = river_excel_file.read_spreadsheet_range(data_t,x_range,y_range)
            for y in self.river_years[:]:
                river_sheet_collected_data[str(y)] = ry[:,count].copy()
                logging.debug(str(y))
            self.river_collected_data[data_t] =  river_sheet_collected_data.copy()    
        logging.debug("End river data collection")
              
        runoff_excel_file = xlsobj.xlsx(self.path_runoff)
        self.runoff_montly_mod = river_excel_file.read_spreadsheet_allrow("monthly",range_montly)
        self.runoff_coordr = river_excel_file.read_spreadsheet_allrow("monthly",range_coord)
        self.nrunoff = len(self.river_coordr[:])
        self.runoff_collected_data = {}
        for data_t in self._mesh_father.input_data.river_data_sheet:
            runoff_sheet_collected_data = {}
            x_range_coord  = [1]
            y_range = range(7,48)
            self.runoff_years = river_excel_file.read_spreadsheet_range(data_t,x_range_coord,y_range,"i")
            count = 0
            x_range = range(2,39)
            ry = runoff_excel_file.read_spreadsheet_range(data_t,x_range,y_range)
            for y in self.river_years[:]:
                runoff_sheet_collected_data[str(y)] = ry[:,count].copy()
            self.runoff_collected_data[data_t] = runoff_sheet_collected_data.copy()
        logging.debug("End runoff data collection")
            
    def _coast_line_mask(self, mask):
        [rows, cols] = np.nonzero(mask)
        # [rows, cols] = find(mask);
        nSea = len(rows);
        coast=mask.astype(np.bool)
        
        for i in range(0,nSea):
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
        x,y = self._mesh_father.tmask[0,1][:].shape
        a = self._coast_line_mask(mask1)
        b = self._coast_line_mask(mask2)
        coast = np.zeros((x,y),dtype=np.bool)
        for i in range(x):
            for j in range(y):
                coast[i,j] = a[i,j] and b[i,j]
        #coast = (self._coast_line_mask(mask1) and self._coast_line_mask(mask2))
        
        loncm = self._mesh_father.nav_lon[coast]
        latcm = self._mesh_father.nav_lat[coast]
        coastline_row_index,coastline_col_index = np.nonzero(coast)
#         coastline_row_index = []
#         coastline_col_index = []
#         for i in coastline_idx: 
#             coastline_row_index.append(i[0])
#             coastline_col_index.append(i[1])
        
        georef4 = np.matrix((coastline_row_index, coastline_col_index, loncm, latcm)).T 
        
        data_types = self._mesh_father.input_data.river_data_sheet
        #n_data_types = len(data_types)
        
        ### river contributes
        georef = np.zeros((self.nrivers,5))
        for jr in range (0,self.nrivers):
            lon_river = self.river_coordr[jr,0]
            lat_river = self.river_coordr[jr,1]
            dist = (loncm-lon_river)**2 + (latcm-lat_river)**2
            ind = np.argmin(dist)
            # w = np.min(dist)
#             print("ind ",ind)
#             print("jr ",jr)
#             print("georef",georef.shape)
            georef[jr,0]=jr
            for i in range(1,5):
                georef[jr,i]=georef4[ind,i-1]
#             print("georef4",np.size(georef4))
        self.river_georef = georef
        
        m=np.zeros((self.nrivers,12))
        self.river_data={}
        
        for data_type in data_types :
            years_data={}
            for ic in self.river_years :
                for r in range(0,self.nrivers-2) :
#                     print(r)
#                     print(self.river_collected_data[data_type][str(ic)].shape)
                    ry = self.river_collected_data[data_type][str(ic)][r]
                    m[r,:] =  (self.river_montly_mod[r,:]/100)*12*ry
                years_data[str(ic)]=m.copy()
            self.river_data[data_type]=years_data.copy() 
        
        ### runoff contributes 
        self.n_coast_cells = len(loncm)
        
        indexes = np.zeros(self.n_coast_cells)
        georef = np.zeros((self.n_coast_cells,5))
        
        for i in range(0,self.n_coast_cells):
            lon_coast_cell = loncm[i]
            lat_coast_cell = latcm[i]
            dist =( (self.runoff_coordr[:,0]-lon_coast_cell)**2
                     + (self.runoff_coordr[:,1]-lat_coast_cell)**2 )
            ind = np.argmin(dist)
            # w = np.min(dist)
            indexes[i] = ind
            georef[i,0]=i
            for ii in range(1,5):
                georef[i,ii]=georef4[i,ii-1]
            #georef[i,:]=np.array(ind ,georef4[i,:])
        self.runoff_georef = georef
        
        
        m=np.zeros((self.nrunoff,12))
        self.runoff_data={}
        
        for data_type in data_types :
            years_data={}
            for ic in self.river_years :
                years_data[str(ic)]=np.zeros((self.n_coast_cells,12)).copy()
                for r in range (0,self.nrunoff-2) :
                    ry = self.runoff_collected_data[data_type][str(ic)][r]
                    m[r,:] =  (self.runoff_montly_mod[r,:]/100)*12*ry
                    ii=indexes==r
                    count = ii.sum()
                    if count > 0:    
                        years_data[str(ic)][ii,:]= npmat.repmat(m[r,:]/count, count, 1)
                       
                    
            self.runoff_data[data_type]=years_data.copy()                    
        
        #sum contributes
    
        for k in range(0,np.size(self.river_georef[0,:])):
            im = self.river_georef[k,2]
            jm = self.river_georef[k,3]
            for i in range(0,self.nrivers):
                if (self.river_georef[i,2] == im and self.river_georef[i,3] == jm) :
                    ii = self.river_georef[i,2]
            for dt in data_types :
                for yr in self.river_years :
                    self.runoff_data[dt][str(yr)][ii,:] = (
                        self.runoff_data[dt][str(yr)][ii,:] +
                        self.river_data[dt][str(yr)][k,:] )
                 
            
        self.river_runoff_data = self.runoff_data


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
        self.gibilterra = lateral_bc(self,self.input_data.file_nutrients)
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
        
        jpk = self.tmask_dimension[1]
        jpj = self.tmask_dimension[2]
        jpi = self.tmask_dimension[3]
        n_coast_cell = self.river.river_georef.shape[0]
        area = self.e1t * self.e2t
        area = area[0,0][:]
        print(area.shape)
        jpt_gib = 4
        nvar_gib = 6
        index = self.bounmesh.idx
        
        for yr in self.river.river_years:
            
            for time in range(1,4):
                name_file = str(yr)+self.gibilterra.season[time]
                for jn in range(self.bounmesh.nudg):
                    aux = self.bounmesh.resto[jn][:]
                    isNudg = np.zeros(aux.shape,dtype=int)
                    
                    for k in range(aux.shape[0]):
                        for j in range(aux.shape[1]):
                            for i in range(aux.shape[2]):
                                if (aux[k,j,i]!=0 and aux[k,j,i] < 1e+19):   # da controllare
                                    isNudg[k,j,i] = 1
                    count = 0
                    npi=np.sum(isNudg);
                    idx=np.zeros(npi); 
                    data=np.zeros(npi);
                    
                    if jn == 1:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.phos[time,jk,jj,ji];
                                        count =  count +1;
                                        
                    if jn == 2:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.ntra[time,jk,jj,ji];
                                        count =  count +1; 
                    
                    if jn == 3:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.dox[time,jk,jj,ji];
                                        count =  count +1;
                    
                    if jn == 4:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.sica[time,jk,jj,ji];
                                        count =  count +1;
                    
                    if jn == 5:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.dic[time,jk,jj,ji];
                                        count =  count +1;
                    if jn == 6:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = self.gibilterra.alk[time,jk,jj,ji];
                                        count =  count +1;
                    if jn == 7:
                        for jk in range(jpk):
                            for jj in range(jpj):
                                for ji in range(jpi):
                                    if isNudg[jk,jj,ji]:
                                        idx[count] = index[jk,jj,ji];
                                        data[count] = 0.0025;
                                        count =  count +1;  

         #  field_name = ['gib_idxt_' vnudg(jn,:) ];
         #  field_data = ['gib_' vnudg(jn,:) ];
         #  S.DIMS.(field_name)=count;
         #  S.(field_name).value=idx;S.(field_name).type='INT' ;
         #  S.(field_data).value=data;S.(field_data).type='DOUBLE' ;
#        filename = ['OUTPUT_DATA/GIB_' timeSTR(I,:) '.nc'];
#     ncwrite(filename,S)

            jpt_riv = 12
            index_riv_a = np.zeros(( 2 * n_coast_cell, 1), dtype = np.int);
            phos_riv_a  = np.zeros(( 2 * n_coast_cell,  jpt_riv));
            ntra_riv_a  = np.zeros(( 2 * n_coast_cell,  jpt_riv));
            sili_riv_a  = np.zeros(( 2 * n_coast_cell,  jpt_riv));
            alka_riv_a  = np.zeros(( 2 * n_coast_cell,  jpt_riv));
            dicc_riv_a  = np.zeros(( 2 * n_coast_cell,  jpt_riv));
            w= 1.0e+12;
            t = 1/(365 * 86400)
            n = 1/14;  
            totN = 0;
            p = 1/31;
            totP = 0;
            s = 1/28;
            totS = 0; 
            totA = 0; 
            totD = 0;
            
            for jc in range(n_coast_cell):
                jc2  = jc + n_coast_cell;
                
                print(self.river.river_georef.shape)
                ji  = self.river.river_georef[jc,2];
                jj  = self.river.river_georef[jc,3];
                Vol2cells = area[jj,ji]*(self.e3t[0,0][:]+self.e3t[0,1][:]);
                cn = w*t*n/Vol2cells;
                cp = w*t*p/Vol2cells;
                cs = w*t*s/Vol2cells; 
                ca = w*t  /Vol2cells; 
                cc = w*t  /Vol2cells; 
                totN = totN + sum(self.river.river_runoff_data["no3_kt_yr"][str(yr)][jc,:],2)/12;
                totP = totP + sum(self.river.river_runoff_data["po4_kt_yr"][str(yr)][jc,:],2)/12;
                totS = totS + sum(self.river.river_runoff_data["sic_kt_yr"][str(yr)][jc,:],2)/12;
                totA = totA + sum(self.river.river_runoff_data["alk_Gmol_yr"][str(yr)][jc,:],2)/12;
                totD = totD + sum(self.river.river_runoff_data["dic_kt_yr"][str(yr)][jc,:],2)/12;
                index_riv_a[jc]  = index[1,jj,ji];
                ntra_riv_a[jc,:] = self.river.river_runoff_data["no3_kt_yr"][str(yr)][jc,:]*cn;
                phos_riv_a[jc,:] = self.river.river_runoff_data["po4_kt_yr"][str(yr)][jc,:]*cp;
                sili_riv_a[jc,:] = self.river.river_runoff_data["sic_kt_yr"][str(yr)][jc,:]*cs;
                alka_riv_a[jc,:] = self.river.river_runoff_data["alk_Gmol_yr"][str(yr)][jc,:]*ca;
                dicc_riv_a[jc,:] = self.river.river_runoff_data["dic_kt_yr"][str(yr)][jc,:]*cc;
                index_riv_a[jc2]  = index[2,jj,ji];
                ntra_riv_a[jc2,:] = self.river.river_runoff_data["no3_kt_yr"][str(yr)][jc,:]*cn;
                phos_riv_a[jc2,:] = self.river.river_runoff_data["po4_kt_yr"][str(yr)][jc,:]*cp;
                sili_riv_a[jc2,:] = self.river.river_runoff_data["sic_kt_yr"][str(yr)][jc,:]*cs;
                alka_riv_a[jc2,:] = self.river.river_runoff_data["alk_Gmol_yr"][str(yr)][jc,:]*ca;
                dicc_riv_a[jc2,:] = self.river.river_runoff_data["dic_kt_yr"][str(yr)][jc,:]*cc;
            print(index)
            print("------------------------------------------------")
            print(index_riv_a)
            ix = np.sort(index_riv_a,1);
            n3n_riv = ntra_riv_a[ix,:];
            n1p_riv = phos_riv_a[ix,:];  
            o3h_riv = alka_riv_a[ix,:];
            o3c_riv = dicc_riv_a[ix,:];
            n5s_riv = sili_riv_a[ix,:]; 
            count_riv = 2 * n_coast_cell;
 
# for I=1:12
#     clear S
#     month= num2str(I,'%02d');
#     timestr = ['2000' month '15-00:00:00'] ;
# 
#     S.DIMS.riv_idxt = count_riv ;
#     S.Attributes    = Attributes;
#     S.Attributes.Type                     = 'Terrestrial Inputs' ;
#     S.Attributes.Month                    = datestr(datenum(timestr,'yyyymmdd-hh:MM:ss'),'mmm') ;
#     S.Attributes.RIV_P_MassBalance_kTON_y = totPriv ;
#     S.Attributes.RIV_N_MassBalance_kTON_y = totNriv ;
#     S.Attributes.RIV_S_MassBalance_kTON_y = totSriv ;
#     S.Attributes.RIV_A_MassBalance_kTON_y = totAriv ;
#     S.Attributes.RIV_D_MassBalance_kTON_y = totDriv ;
#     
#     
#     S.riv_idxt.value = double(idxt_riv); S.riv_idxt.type = 'INT' ; 
#     S.riv_N1p.value  = N1p_riv(:,I)    ; S.riv_N1p.type  = 'DOUBLE' ;
#     S.riv_N3n.value  = N3n_riv(:,I)    ; S.riv_N3n.type  = 'DOUBLE' ;
#     S.riv_N5s.value  = N5s_riv(:,I)    ; S.riv_N5s.type  = 'DOUBLE' ;
#     S.riv_O3c.value  = O3c_riv(:,I)    ; S.riv_O3c.type  = 'DOUBLE' ;
#     S.riv_O3h.value  = O3h_riv(:,I)    ; S.riv_O3h.type  = 'DOUBLE' ;
#     
#     timestr(1:4)=yyyy;
#     filename = ['OUTPUT_DATA/TIN_' timestr '.nc']; 
#     ncwrite(filename,S)
#     

            jpt_atm = 1;
            atm = self.submesh.atm
            count_atm = atm.shape[0];
            index_atm_a = np.zeros(count_atm, dtype = np.int);
            phos_atm_a  = np.zeros((count_atm,jpt_atm));
            ntra_atm_a  = np.zeros((count_atm,jpt_atm));
            w = 1.0e+09;
            t = 1/(365 * 86400);
            totP = 0;
            totN = 0;
            totN_KTy =0; 
            totP_KTy =0; 
            for jn in range(1,count_atm):
                ji  = atm[jn,2];
                jj  = atm[jn,3];
                jk  = 1;
                VolCell1=area[jj,ji]*self.e3t[1];
                totN = totN + atm[jn,6]*VolCell1;
                totP = totP + atm[jn,7]*VolCell1;
                  
                totN_KTy = totN_KTy + atm[jn,6]*(1.e-3/n)*VolCell1;
                totP_KTy = totP_KTy + atm[jn,7]*(1.e-3/p)*VolCell1;
                cn = w*t;
                cp = w*t;
                index_atm_a[jn] = self.bounmesh.idx[jk,jj,ji];
                ntra_atm_a[jn,:] = atm[jn,6]*cn;
                phos_atm_a[jn,:] = atm[jn,7]*cp;

# clear S
# S.DIMS.atm_idxt = count_atm; 
# S.Attributes                          = Attributes; 
# S.Attributes.ATM_P_MassBalance_kTON_y = totP_KTy ; 
# S.Attributes.ATM_N_MassBalance_kTON_y = totN_KTy ; 
# S.Attributes.ATM_P_MassBalance_Mmol_y = totP ; 
# S.Attributes.ATM_N_MassBalance_Mmol_y = totN ;
# S.Attributes.Type                     = 'Atmospherical Inputs' ;
# S.Attributes.Time                     = 'Year' ;
# 
# S.atm_idxt.value = double(index_atm_a ) ; S.atm_idxt.type = 'INT';  
# S.atm_N1p.value  = phos_atm_a           ; S.atm_N1p.type = 'DOUBLE' ;
# S.atm_N3n.value  = ntra_atm_a           ; S.atm_N3n.type = 'DOUBLE' ;
# %S.atm_O2o.value  =  dox_atm_a           ; S.atm_O2o.type = 'DOUBLE' ;
# %S.atm_N5s.value  = sica_atm_a           ; S.atm_N5s.type = 'DOUBLE' ;
# 
# filename = ['OUTPUT_DATA/ATM_' yyyy '0630-00:00:00.nc'];
# ncwrite(filename,S); 

        
        
        
        
        
        
        
        
