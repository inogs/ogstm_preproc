# -*- coding: utf-8 -*-
import numpy as np
import numpy.matlib as npmat
import netCDF4 as nc
from bclib.io_lib import excel_obj as xlsobj
import logging
import code
import matplotlib.pyplot as plt

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
        logging.debug("--Start river data collection")
        work_sheet  = self._mesh_father.input_data.river_data_sheet

        range_montly = range(10,21)
        range_coord  = [2,3]
        force_coord  = [3,4]
        river_excel_file = xlsobj.xlsx(self.path_river)
        river_spreadsheet = {}
        #read xlsx
        river_spreadsheet["monthly"] =  river_excel_file.read_spreadsheet_all("monthly")
        for data_t in self._mesh_father.input_data.river_data_sheet:
            river_spreadsheet[data_t] =  river_excel_file.read_spreadsheet_all(data_t)
        
        #extract data
        self.river_coordr = river_spreadsheet["monthly"][1:,1:3].astype(np.float32)
        self.force_coordr = river_spreadsheet["monthly"][1:,3:5].astype(np.float32)
        self.river_montly_mod = river_spreadsheet["monthly"][1:,9:21].astype(np.float32)
      
        self.nrivers = len(self.river_coordr[:])
        self.river_collected_data = {}
        for data_t in self._mesh_father.input_data.river_data_sheet:
            river_sheet_collected_data = {}
            #self.river_years = river_excel_file.read_spreadsheet_range(data_t,x_range_coord,y_range,"i")
            self.river_years = river_spreadsheet[data_t][0][9:]
            count = 0
            x_range = range(2,41)
            ry = river_spreadsheet[data_t][0][9:]
            count = 8 
            #print(river_spreadsheet[data_t].shape)
            for y in self.river_years[:]:
                count = count + 1
                #print(y,count)
                river_sheet_collected_data[str(y)] =  river_spreadsheet[data_t][1:,count].copy()
            self.river_collected_data[data_t] =  river_sheet_collected_data.copy()
        logging.debug("--End river data collection")

        

        # runoff_excel_file = xlsobj.xlsx(self.path_runoff)
        # self.runoff_montly_mod = river_excel_file.read_spreadsheet_allrow("monthly",range_montly)
        # self.runoff_coordr = river_excel_file.read_spreadsheet_allrow("monthly",range_coord)
        # self.nrunoff = len(self.river_coordr[:])
        # self.runoff_collected_data = {}
        # logging.debug("--Start runoff data collection")
        # for data_t in self._mesh_father.input_data.river_data_sheet:
        #     runoff_sheet_collected_data = {}
        #     x_range_coord  = [1]
        #     y_range = range(7,48)
        #     self.runoff_years = river_excel_file.read_spreadsheet_range(data_t,x_range_coord,y_range,"i")
        #     count = 0
        #     x_range = range(2,39)
        #     ry = runoff_excel_file.read_spreadsheet_range(data_t,x_range,y_range)
        #     for y in self.river_years[:]:
        #         runoff_sheet_collected_data[str(y)] = ry[:,count].copy()
        #     self.runoff_collected_data[data_t] = runoff_sheet_collected_data.copy()
        # logging.debug("--End runoff data collection")


    def _coast_line_mask(self, mask):

        x,y = mask.shape
        coast=np.zeros((x,y),dtype=bool)
        for r in range(1,x-1):
            for c in range(1,y-1):
                expr = (mask[r,c-1] == 0 or mask[r,c+1] == 0 or mask[r-1,c] == 0 or mask[r+1,c] == 0)
                expr = expr or (mask[r-1,c-1] == 0 or mask[r+1,c+1] == 0 or mask[r-1,c+1] == 0 or mask[r+1,c-1] == 0)
                if mask[r,c] == 1 and expr :
                    coast[r,c] = True
        
        return coast

    def map_contribute_on_sea(self):
        """
            Map contributes of terrestrial input on sea
            author : gcoidessa,icelic
        """
        logging.info("Start calc river and runoff calculation")
        mask1 = self._mesh_father.tmask[0,0,:,:]
        coast = self._coast_line_mask(mask1)

        #mask2 = self._mesh_father.tmask[0,1,:,:]
        #x,y = self._mesh_father.tmask[0,0][:].shape
        #b = self._coast_line_mask(mask2)
        #coast = np.zeros((x,y),dtype=np.bool)
        #for i in range(x):
        #    for j in range(y):
        #        coast[i,j] = a[i,j] and b[i,j]

        loncm = self._mesh_father.nav_lon[coast]
        latcm = self._mesh_father.nav_lat[coast]
        coastline_row_index,coastline_col_index = np.nonzero(coast)

        georef4 = np.matrix((coastline_row_index, coastline_col_index, loncm, latcm)).T

        data_types = self._mesh_father.input_data.river_data_sheet
       
        #code.interact(local=locals())

        ### river contributes
        georef = np.zeros((self.nrivers,5))
        for jr in range (self.nrivers):
            #print(jr)
            lon_river = self.river_coordr[jr,0]
            lat_river = self.river_coordr[jr,1]
            #print(loncm,"-",lon_river,"+",latcm,"-",lat_river)
            dist = (loncm-lon_river)**2 + (latcm-lat_river)**2
            #print("dist=",dist, dist.shape)
            ind = np.argmin(dist)
            #print("ind =",ind,georef4.shape)
            #print(georef4[ind,:])
            georef[jr,0]=jr
            for i in range(1,5):
                georef[jr,i]=georef4[ind,i-1]
            #print(self.force_coordr[jr,:])
            #force cordinates
            if self.force_coordr[jr,0] != -1 and self.force_coordr[jr,1] != -1:
                georef[jr,1]=self.force_coordr[jr,1]
                georef[jr,2]=self.force_coordr[jr,0]
            #print(georef[jr,:])
        self.river_georef = georef
        
        m=np.zeros((self.nrivers,12))
        self.river_data={}

        for data_type in data_types :
            years_data={}
            for ic in self.river_years :
                for r in range(0,self.nrivers-2):
                    ry = self.river_collected_data[data_type][str(ic)][r]
                    m[r,:] =  (self.river_montly_mod[r,:]/100)*12*ry
                years_data[str(ic)]=m.copy()
            self.river_data[data_type]=years_data.copy()

        ### runoff contributes
        # self.n_coast_cells = len(loncm)

        # indexes = np.zeros(self.n_coast_cells)
        # georef = np.zeros((self.n_coast_cells,5))

        # for i in range(0,self.n_coast_cells):
        #     lon_coast_cell = loncm[i]
        #     lat_coast_cell = latcm[i]
        #     dist =( (self.runoff_coordr[:,0]-lon_coast_cell)**2
        #              + (self.runoff_coordr[:,1]-lat_coast_cell)**2 )
        #     ind = np.argmin(dist)
        #     indexes[i] = ind
        #     georef[i,0]=i
        #     for ii in range(1,5):
        #         georef[i,ii]=georef4[i,ii-1]
        # self.runoff_georef = georef


        # m=np.zeros((self.nrunoff,12))
        # self.runoff_data={}

        # for data_type in data_types :
        #     years_data={}
        #     for ic in self.river_years :
        #         years_data[str(ic)]=np.zeros((self.n_coast_cells,12)).copy()
        #         for r in range (0,self.nrunoff-2) :
        #             ry = self.runoff_collected_data[data_type][str(ic)][r]
        #             m[r,:] =  (self.runoff_montly_mod[r,:]/100)*12*ry
        #             ii=indexes==r
        #             count = ii.sum()
        #             if count > 0:
        #                 years_data[str(ic)][ii,:]= npmat.repmat(m[r,:]/count, count, 1)


        #     self.runoff_data[data_type]=years_data.copy()

        # #sum contributes

        # for k in range(0,np.size(self.river_georef[0,:])):
        #     im = self.river_georef[k,2]
        #     jm = self.river_georef[k,3]
        #     for i in range(0,self.nrivers):
        #         if (self.river_georef[i,2] == im and self.river_georef[i,3] == jm) :
        #             ii = self.river_georef[i,2]
        #     for dt in data_types :
        #         for yr in self.river_years :
        #             self.runoff_data[dt][str(yr)][ii,:] = (
        #                 self.runoff_data[dt][str(yr)][ii,:] +
        #                 self.river_data[dt][str(yr)][k,:] )


        # self.river_runoff_data = self.runoff_data
        logging.info("End calc river and runoff calculation")
