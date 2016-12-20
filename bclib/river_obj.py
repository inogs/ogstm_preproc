# -*- coding: utf-8 -*-
import numpy as np
from bclib import excel_obj as xlsobj
import logging
import netCDF4


class river_data:

    def __init__(self,conf):
        '''
        Reads data from excel in 
        river_collected_data, a dict of dicts
        '''
        self.georef = None
        logging.info("--Start river data collection")
        sheet_list=conf.river_data_sheet
        river_excel_file = xlsobj.xlsx(conf.file_river)
        river_spreadsheet = {}
        #read xlsx into a dict
        river_spreadsheet["monthly"] =  river_excel_file.read_spreadsheet_all("monthly")
        for sheet in sheet_list:
            river_spreadsheet[sheet] =  river_excel_file.read_spreadsheet_all(sheet)

        #extract data
        self.river_coordr = river_spreadsheet["monthly"][1:,1:3].astype(np.float32)
        self.force_coordr = river_spreadsheet["monthly"][1:,3:5].astype(np.int32)
        self.river_montly_mod = river_spreadsheet["monthly"][1:,9:21].astype(np.float32)

        self.nrivers = len(self.river_coordr[:])
        self.river_collected_data = {}
        for sheet in sheet_list:
            river_sheet_collected_data = {}
            self.river_years = river_spreadsheet[sheet][0][9:]
            for iyear, y in enumerate(self.river_years):
                river_sheet_collected_data[str(y)] =  river_spreadsheet[sheet][1:,iyear+9].copy()
            self.river_collected_data[sheet] =  river_sheet_collected_data.copy()
        logging.info("--End river data collection")



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
        logging.info("Start river position calculation")
        there_are_free_points=np.any(self.force_coordr==-1)
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
            if self.force_coordr[jr,0] != -1 and self.force_coordr[jr,1] != -1:
                georef['indLon'][jr]=self.force_coordr[jr,0]
                georef['indLat'][jr]=self.force_coordr[jr,1]
                #if(mask1[georef[jr,1]-1,georef[jr,2]-1] == 0):
                #    print("RIVER ON THE LAND")
                #    print(georef[jr,:])
            else: # TO BE TESTED
                lon_river = self.river_coordr[jr,0]
                lat_river = self.river_coordr[jr,1]
                dist = (loncm-lon_river)**2 + (latcm-lat_river)**2
                ind = np.argmin(dist)
                georef[jr,0]=jr
                for i in range(1,5): georef[jr,i]=georef4[ind,i-1]
                georef[jr,1]=self.georef[jr,1]+1
                georef[jr,2]=self.georef[jr,2]+1

        self.georef = georef
        logging.info("End river position calculation")


    def modularize(self,conf):

        '''
         Applyies monthly modulation to yearly data
         Generates a new field river_data, 
         it is a dict of dicts
         river_data={'sheet_name':{'2009': (nRivers,12) np.array } ; ... }
        '''

        logging.info("Start river Modularization")
        m=np.zeros((self.nrivers,12))
        river_data={}
        sheet_list = conf.river_data_sheet
        for sheet in sheet_list :
            years_data={}
            for ic in self.river_years :
                for r in range(0,self.nrivers-2):
                    ry = self.river_collected_data[sheet][str(ic)][r]
                    m[r,:] =  (self.river_montly_mod[r,:]/100)*12*ry
                years_data[str(ic)]=m.copy()
            river_data[sheet]=years_data.copy()

        self.river_data = river_data
        logging.info("End river Modularization")

    def gen_boun_indexes(self,boun_indexes):
        '''
        boun_indexes is the index array of bounmask.nc
        Bounmask indexes are supposed start from one --> the index_riv_a array is good for fortran
        excel data i,j  are indexes starting from one
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
        
        return N,P,S,A,D
    
    
    
    def conversion(self,N,P,S,A,D):
        '''
        Returns
        * N * in mmol/s
        * P * in mmol/s
        * S * in mmol/s
        * A * in   mg/s
        * D * in   mg/s
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
        return N*cn, P*cp, S*cs,A*ca, D*cc

    def generate_monthly_files(self,conf,mask, idxt_riv, positions):
        '''
        Generates TIN files for every year and every month.

        positions is the same of georef, then it starts from zero
        '''
        logging.info("Start non climatological TIN file generation ")
        Area=np.zeros((self.nrivers,),np.float)

        for jr in range(self.nrivers):
            ji = self.georef['indLon'][jr]-1
            jj = self.georef['indLat'][jr]-1
            Area[jr] = mask.area[jj,ji]

        start_year=conf.simulation_start_time
        end___year=conf.simulation_end_time
        for year in range(start_year,end___year):
            for month in range(1,13):
                filename = conf.dir_out+"/TIN_%d%02d15-00:00:00.nc" %(year, month)
                N,P,S,A,D = self.get_monthly_data(str(year), month)
                N,P,S,A,D = self.conversion(N, P, S, A, D)
                self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, idxt_riv, positions)
        logging.info("End non climatological TIN file generation")

    def generate_climatological_monthly_files(self,conf,mask,idxt_riv, positions):
        '''
        Generates 12 TIN_yyyy*nc files
        '''
        logging.info("Start climatological TIN file generation ")
        Area=np.zeros((self.nrivers,),np.float)

        for jr in range(self.nrivers):
            ji = self.georef['indLon'][jr]-1
            jj = self.georef['indLat'][jr]-1
            Area[jr] = mask.area[jj,ji]

        year="yyyy"
        for month in range(1,13):
            filename = conf.dir_out+"/TIN_yyyy%02d15-00:00:00.nc" %(month)
            N,P,S,A,D = self.get_monthly_data(str(year), month)
            N,P,S,A,D = self.conversion(N, P, S, A, D)
            self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, idxt_riv, positions)
        logging.info("Start climatological TIN file generation ")
                
                #N/A, P/A
    def dump_file(self,filename,N,P,S,A,D,idxt_riv,positions):
        ncfile = netCDF4.Dataset(filename, 'w')
        ncfile.createDimension("riv_idxt",self.nrivers)
        ncfile.createDimension("cords",3)
        riv_idxt_riv = ncfile.createVariable('riv_idxt', 'i4', ('riv_idxt',))
        riv_pos_riv = ncfile.createVariable('position', 'i4', ('riv_idxt','cords'))
        riv_a_n3n = ncfile.createVariable('riv_N3n', 'f4', ('riv_idxt',))
        riv_a_n1p = ncfile.createVariable('riv_N1p', 'f4', ('riv_idxt',))
        riv_a_n5s = ncfile.createVariable('riv_N5s', 'f4', ('riv_idxt',))
        riv_a_o3c = ncfile.createVariable('riv_O3c', 'f4', ('riv_idxt',))
        riv_a_o3h = ncfile.createVariable('riv_O3h', 'f4', ('riv_idxt',))
        
        riv_idxt_riv[:] = idxt_riv[:]
        riv_pos_riv[:,:] = positions
        riv_a_n3n[:] = N
        riv_a_n1p[:] = P
        riv_a_n5s[:] = S
        riv_a_o3c[:] = D
        riv_a_o3h[:] = A
        ncfile.close()
                