# -*- coding: utf-8 -*-
import numpy as np
from bclib import excel_reader
import logging
import netCDF4


class river():

    def __init__(self,conf):
        '''
        Reads data from excel in 
        xls_data, a dict of dicts

        Example
        import config as conf
        R=river(conf)
        a = R.xls_data['DIP_KTperYR_NOBLS']['2004']
        '''
        self.georef = None
        logging.info("Reading from xls file...")
        sheet_list=conf.river_data_sheet
        river_excel_file = excel_reader.xlsx(conf.file_river)
        river_spreadsheet = {}
        #read xlsx into a dict
        river_spreadsheet["monthly"] =  river_excel_file.read_spreadsheet_all("monthly")
        for sheet in sheet_list:
            river_spreadsheet[sheet] =  river_excel_file.read_spreadsheet_all(sheet)

        #extract data
        self.lon = river_spreadsheet["monthly"][1:,1].astype(np.float32)
        self.lat = river_spreadsheet["monthly"][1:,2].astype(np.float32)
        self.forced_coord = river_spreadsheet["monthly"][1:,3:5].astype(np.int32)
        self.monthly_mod  = river_spreadsheet["monthly"][1:,9:21].astype(np.float32)

        self.nrivers = len(self.lon)
        self.xls_data = {}
        for sheet in sheet_list:
            river_sheet_collected_data = {}
            self.available_years = river_spreadsheet[sheet][0][9:]
            for iyear, y in enumerate(self.available_years):
                river_sheet_collected_data[str(y)] =  river_spreadsheet[sheet][1:,iyear+9].copy()
            self.xls_data[sheet] =  river_sheet_collected_data.copy()
        logging.info("Done")



    def _coast_line_mask(self, mask):
        '''
        mask is a 2d array of booleans
        Returns:
        coast:  a 2d array of booleans
        '''

        x,y = mask.shape
        coast=np.zeros((x,y),dtype=bool)
        for r in range(1,x-1):
            for c in range(1,y-1):
                expr = (mask[r,c-1] == 0 or mask[r,c+1] == 0 or mask[r-1,c] == 0 or mask[r+1,c] == 0)
                expr = expr or (mask[r-1,c-1] == 0 or mask[r+1,c+1] == 0 or mask[r-1,c+1] == 0 or mask[r+1,c-1] == 0)
                if mask[r,c] == 1 and expr :
                    coast[r,c] = True

        return coast

    def gen_map_indexes(self,mask):
        """
        Generates the "georef" field
          a (nRivers,) numpy array
        indLon, indLat, lonmesh, latmesh

        Warning:
        * for forced coordinates indLon and indLat are starting from one
        * the others are not tested
        """
        logging.info("River position calculation: start ")
        there_are_free_points=np.any(self.forced_coord==-1)
        if there_are_free_points:
            mask1 = mask.mask_at_level(0)
            coast = self._coast_line_mask(mask1)
    
            loncm = mask.xlevels[coast]
            latcm = mask.ylevel[coast]
            coastline_row_index,coastline_col_index = np.nonzero(coast)
    
            georef4 = np.matrix((coastline_row_index, coastline_col_index, loncm, latcm)).T


        ### river contributes
        georef = np.zeros((self.nrivers,),dtype=[('indLon',np.int),('indLat',np.int),('lonmesh',np.float32),('latmesh',np.float32)])
        for jr in range (self.nrivers):
            if self.forced_coord[jr,0] != -1 and self.forced_coord[jr,1] != -1:
                georef['indLon'][jr]=self.forced_coord[jr,0]
                georef['indLat'][jr]=self.forced_coord[jr,1]
                #if(mask1[georef[jr,1]-1,georef[jr,2]-1] == 0):
                #    print("RIVER ON THE LAND")
                #    print(georef[jr,:])
            else: # TO BE TESTED
                dist = (loncm-self.lon[jr])**2 + (latcm-self.lat[jr])**2
                ind = np.argmin(dist)
                georef[jr,0]=jr
                for i in range(1,5): georef[jr,i]=georef4[ind,i-1]
                georef[jr,1]=self.georef[jr,1]+1
                georef[jr,2]=self.georef[jr,2]+1

        self.georef = georef
        logging.info("River position calculation: done ")


    def modularize(self,conf):

        '''
         Applyies monthly modulation to yearly data
         Generates a new field river_data,
         it is a dict of dicts
         river_data={'sheet_name':{'2009': (nRivers,12) np.array } ; ... }
         Modularize factors are read in "monthly" sheet,
         they are pure factors, not percent.
         Modularization takes in account weight of eanch month.
        '''

        logging.info("River Modularization")
        days_in_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        m=np.zeros((self.nrivers,12),np.float32)
        river_data={}
        sheet_list = conf.river_data_sheet
        for sheet in sheet_list :
            years_data={}
            for year in self.available_years:
                yearstr=str(year)
                for r in range(self.nrivers):
                    ry = self.xls_data[sheet][yearstr][r]
                    m[r,:] =  (self.monthly_mod[r,:]*365/days_in_month)*ry
                years_data[yearstr]=m.copy()
            river_data[sheet]=years_data.copy()

        self.river_data = river_data

    def gen_boun_indexes(self,boun_indexes):
        '''
        For each river, generates its bounmask index
        Arguments
        * boun_indexes * is a 3d array of integers, the "index" array of bounmask.nc
        Returns
        * idxt      * integer 1d array (nrivers)   with rivers index in bounmask
                      indexes start from one, is ok for fortran
        * positions * integer 2d array (nrivers,3) with river position jk,jj,ji in meshmask
                      indexes start from zero, is ok for python

        External conditions:
        * bounmask.nc['index'] values are supposed start from one
        * excel data i,j  are indexes starting from one
        '''

        position    = np.zeros((self.nrivers,3), dtype = np.int);
        index_riv_a = np.zeros((self.nrivers,) , dtype = np.int);
        for jr in range(self.nrivers):
            jj  = self.georef[jr]['indLat']-1;
            ji  = self.georef[jr]['indLon']-1;
            index_riv_a[jr]  = boun_indexes[0,jj,ji];
            assert index_riv_a[jr] > 0
            position[jr] = [0,jj,ji]
        idxt_riv = index_riv_a; #no sort at the moment
        return idxt_riv, position

    def total(self,year):
        '''
        Returns:
        * nitrate * kT/yr 
        TOTALMENTE INUTILE ADESSO
        '''
 
        totN = self.river_data["DIN_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        totP = self.river_data["DIP_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        totS = self.river_data["DIS_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        totA = self.river_data["ALK_GmolperYR_NOBLS"][str(year)].sum(axis=1)/12;
        totD = self.river_data["DIC_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        return totN,totP,totS,totA,totD


    def get_monthly_data(self,yearstr,month):
        '''
        Returns monthly data, read from excel and modularized
        Arguments: 
         *yearstr* year as strig
         *month*  integer from 1 to 12
         
         Returns:
         N,P,S,A,D : (nRivers, ) numpy arrays
        
        '''
        N = self.river_data["DIN_KTperYR_NOBLS"  ][yearstr][:,month-1]
        P = self.river_data["DIP_KTperYR_NOBLS"  ][yearstr][:,month-1]
        S = self.river_data["DIS_KTperYR_NOBLS"  ][yearstr][:,month-1]
        A = self.river_data["ALK_GmolperYR_NOBLS"][yearstr][:,month-1]
        D = self.river_data["DIC_KTperYR_NOBLS"  ][yearstr][:,month-1]
        O = self.river_data["O2o_GmolperYR_NOBLS"][yearstr][:,month-1]
        return N,P,S,A,D,O
    
    
    
    def conversion(self,N,P,S,A,D,O):
        '''
        Performs conversion of all variables
        from KT/y to mmol/s or mg/s

        Arguments:
         Nitrate, Phosphate, Silicate, Alcalinity, Dissoved Inorganic Carbon, Dissolved_Oxygen arrays

        Returns
        * N * in mmol/s
        * P * in mmol/s
        * S * in mmol/s
        * A * in   mg/s
        * D * in   mg/s
        * O * in mmol/s
         as (nRivers,) numpy array

        Data are not ready for model, they have to be divided by cell area or cell volume.
        ''' 
        w= 1.0e+12;
        t = 1./(365 * 86400)
        n = 1./14;
        p = 1./31;
        s = 1./28;
        cn = w*t*n
        cp = w*t*p
        cs = w*t*s
        ca = w*t  
        cc = w*t    
        return N*cn, P*cp, S*cs,A*ca, D*cc, O*cc

    def generate_monthly_files(self,conf,mask, idxt_riv, positions):
        '''
        Generates TIN files for every year and every month.

        positions is the same of georef, then it starts from zero
        '''
        logging.info("Non climatological TIN file generation : start ")
        Area=np.zeros((self.nrivers,),np.float)

        for jr in range(self.nrivers):
            ji = self.georef['indLon'][jr]-1
            jj = self.georef['indLat'][jr]-1
            Area[jr] = mask.area[jj,ji]

        start_year=conf.simulation_start_time
        end___year=conf.simulation_end_time
        for year in range(start_year,end___year):
            for month in range(1,13):
                filename = conf.dir_out + "TIN_%d%02d15-00:00:00.nc" %(year, month)
                N,P,S,A,D,O = self.get_monthly_data(str(year), month)
                N,P,S,A,D,O = self.conversion(N, P, S, A, D, O)
                self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, O/Area, idxt_riv, positions)
        logging.info("Non climatological TIN file generation : done")

    def generate_climatological_monthly_files(self,conf,mask,idxt_riv, positions):
        '''
        Generates 12 TIN_yyyy*nc files
        '''
        logging.info("Climatological TIN file generation : start")
        Area=np.zeros((self.nrivers,),np.float)

        for jr in range(self.nrivers):
            ji = self.georef['indLon'][jr]-1
            jj = self.georef['indLat'][jr]-1
            Area[jr] = mask.area[jj,ji]

        year="yyyy"
        for month in range(1,13):
            filename = conf.dir_out+"/TIN_yyyy%02d15-00:00:00.nc" %(month)
            N,P,S,A,D,O = self.get_monthly_data(str(year), month)
            N,P,S,A,D,O = self.conversion(N, P, S, A, D, O)
            self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, O/Area, idxt_riv, positions)
        logging.info("Climatological TIN file generation : done")
                
    def dump_file(self,filename,N,P,S,A,D,O,idxt_riv,positions):
        '''
          Writes the single TIN file
          Variables are dumped as they are, all but positions (incremented by one)
        '''
        ncfile = netCDF4.Dataset(filename, 'w')
        ncfile.createDimension("riv_idxt",self.nrivers)
        ncfile.createDimension("coords",3)
        riv_idxt_riv = ncfile.createVariable('riv_idxt', 'i4', ('riv_idxt',))
        riv_pos   = ncfile.createVariable('position', 'i4', ('riv_idxt','coords'))
        setattr(riv_pos,'order',"k,j,i")
        setattr(ncfile,'Units','mmol or mg /(s*m2) ')
        riv_a_n3n = ncfile.createVariable('riv_N3n', 'f4', ('riv_idxt',))
        riv_a_n1p = ncfile.createVariable('riv_N1p', 'f4', ('riv_idxt',))
        riv_a_n5s = ncfile.createVariable('riv_N5s', 'f4', ('riv_idxt',))
        riv_a_o3c = ncfile.createVariable('riv_O3c', 'f4', ('riv_idxt',))
        riv_a_o3h = ncfile.createVariable('riv_O3h', 'f4', ('riv_idxt',))
        riv_a_O2o = ncfile.createVariable('riv_O2o', 'f4', ('riv_idxt',))
        
        riv_idxt_riv[:] = idxt_riv[:]
        riv_pos[:,:] = positions+1
        riv_a_n3n[:] = N
        riv_a_n1p[:] = P
        riv_a_n5s[:] = S
        riv_a_o3c[:] = D
        riv_a_o3h[:] = A
        riv_a_O2o[:] = O
        ncfile.close()

    def load_from_file(self,filename,var):
        ''' Useful for check/debug '''
        dset = netCDF4.Dataset(filename, 'r')
        M = np.array(dset[var])
        dset.close()
        return M
