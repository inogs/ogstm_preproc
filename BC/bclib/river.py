# -*- coding: utf-8 -*-

# See /gss/gss_work/DRES_OGS_BiGe/gcossarini/River_Generation to know about excel generation

import numpy as np
from bclib import excel_reader
import logging
import netCDF4

def conversion(var):
    '''
    from mmol or mg to KTONS(N1p,N3n,N5s,O3c,O3h)  or Gmol (O2o)
    from nmol to KTONS (Hg, MMHg)
    '''
    Conversion={}
    w= 1.0e-12
    n = 14
    p = 31
    s = 28
    h = 200590000    #/200.59*1000000.
    Conversion['N1p'] =w*p
    Conversion['N3n'] =w*n
    Conversion['N5s'] =w*s
    Conversion['O3h'] =w
    Conversion['O3c'] =w
    Conversion['O2o'] =w
    Conversion['Hg2'] =w*h
    Conversion['MMHg'] =w*h
    return Conversion[var]

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
                river_sheet_collected_data[str(int(y))] =  river_spreadsheet[sheet][1:,iyear+9].astype(np.float64)
            self.xls_data[sheet] =  river_sheet_collected_data.copy()
        self.xls_data['monthly'] = river_spreadsheet["monthly"]
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

    def gen_map_indexes(self,mask, nspread=1, delta=5):
        """
        Arguments:
        * mask    * a mask object
        * nspread * integer, the number of horizontal cells in which a river will be placed
                  * the nearest "nspread" cells to the corresponding river cell will be choosen
        * delta   * integer, the dimension in cells (the half side)
                    of a square box to search for nspread cells.

        Generates the "georef" field
          a (nRivers,) numpy array with indLon, indLat, lonmesh, latmesh
        or the "georef_spread" field
           with the same meaning

        Warning:
        * for forced coordinates indLon and indLat are starting from one
        * the others are not tested
        """
        logging.info("River position calculation: start ")
        self.nspread=nspread
        there_are_free_points=np.any(self.forced_coord==-1)
        if there_are_free_points:
            mask1 = mask.mask_at_level(0)
            coast = self._coast_line_mask(mask1)
    
            loncm = mask.xlevels[coast]
            latcm = mask.ylevels[coast]
            coastline_row_index,coastline_col_index = np.nonzero(coast)


        ### river contributes
        georef = np.zeros((self.nrivers,),                 dtype=[('indLon',int),('indLat',int),('lonmesh',np.float32),('latmesh',np.float32)])
        georef_spread=np.zeros((self.nrivers*self.nspread),dtype=[('indLon',int),('indLat',int),('lonmesh',np.float32),('latmesh',np.float32)])
        for jr in range (self.nrivers):
            if self.forced_coord[jr,0] != -1 and self.forced_coord[jr,1] != -1:
                georef['indLon'][jr]=self.forced_coord[jr,0]
                georef['indLat'][jr]=self.forced_coord[jr,1]
                #if(mask1[georef[jr,1]-1,georef[jr,2]-1] == 0):
                #    print("RIVER ON THE LAND")
                #    print(georef[jr,:])
            else:
                dist = (loncm-self.lon[jr])**2 + (latcm-self.lat[jr])**2
                ind = np.argmin(dist)
                indlon = coastline_col_index[ind]+1
                indlat = coastline_row_index[ind]+1
                if (self.nspread==1):
                    georef['lonmesh'][jr]=loncm[ind]
                    georef['latmesh'][jr]=latcm[ind]

                    georef['indLon'][jr] = coastline_col_index[ind]+1
                    georef['indLat'][jr] = coastline_row_index[ind]+1
                else:
                    lonmesh=loncm[ind]
                    latmesh=latcm[ind]
                    _, jpj, jpi=mask.shape
                    I,J=np.meshgrid(np.arange(jpi), np.arange(jpj))
                    localmask     =mask1[indlat-delta:indlat+delta,indlon-delta:indlon+delta]
                    local_X=mask.xlevels[indlat-delta:indlat+delta,indlon-delta:indlon+delta]
                    local_Y=mask.ylevels[indlat-delta:indlat+delta,indlon-delta:indlon+delta]
                    local_I=           I[indlat-delta:indlat+delta,indlon-delta:indlon+delta]
                    local_J=           J[indlat-delta:indlat+delta,indlon-delta:indlon+delta]
                    nPoints=localmask.sum()
                    NEAREST=np.zeros((nPoints,),dtype=[('I',np.int32), ('J',np.int32),('lon',np.float32), ('lat',np.float32), ('dist',np.float32)])
                    NEAREST['I'   ]=local_I[localmask]
                    NEAREST['J'   ]=local_J[localmask]
                    NEAREST['lon' ]=local_X[localmask]
                    NEAREST['lat' ]=local_Y[localmask]
                    NEAREST['dist'] = (NEAREST['lon']-lonmesh)**2 + (NEAREST['lat']-latmesh)**2
                    NEAREST_ORDERED=np.sort(NEAREST,order='dist')
                    for k_spread in range(self.nspread):
                        ind_spread= jr*self.nspread + k_spread
                        georef_spread['indLon' ][ind_spread]=NEAREST_ORDERED['I'  ][k_spread]+1
                        georef_spread['indLat' ][ind_spread]=NEAREST_ORDERED['J'  ][k_spread]+1
                        georef_spread['lonmesh'][ind_spread]=NEAREST_ORDERED['lon'][k_spread]
                        georef_spread['latmesh'][ind_spread]=NEAREST_ORDERED['lat'][k_spread]
        if (self.nspread==1):
            self.georef = georef
        else:
            self.georef_spread = georef_spread
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
                yearstr=str(int(year))
                for r in range(self.nrivers):
                    ry = self.xls_data[sheet][yearstr][r]
                    m[r,:] =  (self.monthly_mod[r,:]*365/days_in_month)*ry
                years_data[yearstr]=m.copy()
            river_data[sheet]=years_data.copy()

        self.river_data = river_data

    def get_coordinates(self,mask):
        '''
        Arguments:
        * mask    * a mask object
        Returns:
           * LAT *  numpy ndarray
           * LON *  numpy ndarray
        '''
        LAT = np.zeros((self.nrivers), np.float32)
        LON = np.zeros((self.nrivers), np.float32)
        for jr in range(self.nrivers):
            jj  = self.georef[jr]['indLat']-1
            ji  = self.georef[jr]['indLon']-1
            LAT[jr] = mask.ylevels[jj,ji]
            LON[jr] = mask.xlevels[jj,ji]

        return LAT, LON



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

        position    = np.zeros((self.nrivers,3), dtype = int);
        index_riv_a = np.zeros((self.nrivers,) , dtype = int);
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
        totH = self.river_data["HgII_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        totM = self.river_data["MMHg_KTperYR_NOBLS"  ][str(year)].sum(axis=1)/12;
        return totN,totP,totS,totA,totD,totH,totM


    def spread_array(self,array):
        assert len(array)==self.nrivers
        long_array=np.zeros((self.nrivers*self.nspread))
        for jr in range(self.nrivers):
            for k in range(self.nspread):
                ind = jr*self.nspread + k
                long_array[ind] = array[jr]/self.nspread
        return long_array

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
        H = self.river_data["HgII_KTperYR_NOBLS"][yearstr][:,month-1]
        M = self.river_data["MMHg_KTperYR_NOBLS"][yearstr][:,month-1]
        DOC  = self.river_data["DOC_KTperYR_NOBLS"][yearstr][:,month-1]
        CDOM = self.river_data["CDOM_KTperYR_NOBLS"][yearstr][:,month-1]*4.0

        if (self.nspread==1):
            return N,P,S,A,D,O,H,M,DOC,CDOM
        else:
            N_long = self.spread_array(N)
            P_long = self.spread_array(P)
            S_long = self.spread_array(S)
            A_long = self.spread_array(A)
            D_long = self.spread_array(D)
            O_long = self.spread_array(O)
            H_long = self.spread_array(H)
            M_long = self.spread_array(M)
            DOC_long = self.spread_array(DOC)
            CDOM_long = self.spread_array(CDOM)
            return N_long, P_long, S_long, A_long, D_long, O_long,H_long,M_long,DOC_long, CDOM_long


    
    
    def conversion(self,N,P,S,A,D,O,H,M,DOC,CDOM):
        '''
        Performs conversion of all variables
        from KT/y to mmol/s or mg/s or nmol/s

        Arguments:
         Nitrate, Phosphate, Silicate, Alcalinity, Dissoved Inorganic Carbon, Dissolved_Oxygen
         Dissolved Organic Carbon, CDOM
        arrays

        Returns
        * N * in mmol/s
        * P * in mmol/s
        * S * in mmol/s
        * A * in   mg/s
        * D * in   mg/s
        * O * in mmol/s
        * DOC *    mg/s
        * CDOM *   mg/s
        * Hg *   nmol/s
        * MeHg * nmol/s
         as (nRivers,) numpy array

        Data are not ready for model, they have to be divided by cell area or cell volume.
        ''' 
        w= 1.0e+12;
        w2= 1.0e+6;
        t = 1./(365 * 86400)
        n = 1./14;
        p = 1./31;
        s = 1./28;
        h = 1./200.59;
        cn = w*t*n
        cp = w*t*p
        cs = w*t*s
        ca = w*t  
        cc = w*t    
        ch = w*w2*h*t    
        return N*cn, P*cp, S*cs,A*ca, D*cc, O*cc,H*ch,M*ch, DOC*cc, CDOM*cc

    def generate_monthly_files(self,conf,mask):#, idxt_riv, positions):
        '''
        Generates TIN files for every year and every month.

        positions is the same of georef, then it starts from zero
        '''
        logging.info("Non climatological TIN file generation : start ")

        if (self.nspread==1):
            Area=np.zeros((self.nrivers,),float)
            for jr in range(self.nrivers):
                ji = self.georef['indLon'][jr]-1
                jj = self.georef['indLat'][jr]-1
                Area[jr] = mask.area[jj,ji]
        else:
            Area = np.zeros((self.nrivers*self.nspread,),float)
            for jr in range(self.nrivers*self.nspread):
                ji = self.georef_spread['indLon'][jr]-1
                jj = self.georef_spread['indLat'][jr]-1
                Area[jr] = mask.area[jj,ji]

        start_year=conf.simulation_start_time
        end___year=conf.simulation_end_time
        for year in range(start_year,end___year):
            for month in range(1,13):
                filename = conf.dir_out + "TIN_%d%02d15-00:00:00.nc" %(year, month)
                N,P,S,A,D,O,H,M,DOC,CDOM = self.get_monthly_data(str(year), month)
                N,P,S,A,D,O,H,M,DOC, CDOM = self.conversion(N, P, S, A, D, O,H,M,DOC, CDOM)
                self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, O/Area,H/Area, M/Area, DOC/Area, CDOM/Area, mask)
        logging.info("Non climatological TIN file generation : done")

    def generate_climatological_monthly_files(self,conf,mask): #idxt_riv, positions):
        '''
        Generates 12 TIN_yyyy*nc files
        '''
        logging.info("Climatological TIN file generation : start")
        Area=np.zeros((self.nrivers,),np.float32)

        if (self.nspread==1):
            Area=np.zeros((self.nrivers,),np.float32)
            for jr in range(self.nrivers):
                ji = self.georef['indLon'][jr]-1
                jj = self.georef['indLat'][jr]-1
                Area[jr] = mask.area[jj,ji]
        else:
            Area = np.zeros((self.nrivers*self.nspread,),np.float32)
            for jr in range(self.nrivers*self.nspread):
                ji = self.georef_spread['indLon'][jr]-1
                jj = self.georef_spread['indLat'][jr]-1
                Area[jr] = mask.area[jj,ji]

        year="yyyy"
        for month in range(1,13):
            filename = conf.dir_out+"/TIN_yyyy%02d15-00:00:00.nc" %(month)
            N,P,S,A,D,O,H,M,DOC,CDOM = self.get_monthly_data(str(year), month)
            N,P,S,A,D,O,H,M,DOC,CDOM = self.conversion(N, P, S, A, D, O,H,M, DOC,CDOM)
            self.dump_file(filename, N/Area, P/Area, S/Area, A/Area, D/Area, O/Area, H/Area,M/Area,DOC/Area, CDOM/Area, mask)
        logging.info("Climatological TIN file generation : done")
                

    def dump_file_old(self,filename,N,P,S,A,D,O,H,M, idxt_riv,positions):
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
        riv_a_Hg2 = ncfile.createVariable('riv_Hg2', 'f4', ('riv_idxt',))
        riv_a_MMHg = ncfile.createVariable('riv_MMHg', 'f4', ('riv_idxt',))

        riv_idxt_riv[:] = idxt_riv[:]
        riv_pos[:,:] = positions+1
        riv_a_n3n[:] = N
        riv_a_n1p[:] = P
        riv_a_n5s[:] = S
        riv_a_o3c[:] = D
        riv_a_o3h[:] = A
        riv_a_O2o[:] = O
        riv_a_Hg2[:] = H
        riv_a_MMHg[:] = M
        ncfile.close()

    def dump_file(self,filename,N,P,S,A,D,O,H,M,DOC, CDOM, mask):
        '''
          Writes the single TIN file
          Variables are dumped as they are, all but positions (incremented by one)
        '''
        _,jpj, jpi = mask.shape
        ncfile = netCDF4.Dataset(filename, 'w')
        ncfile.createDimension('lon',jpi)
        ncfile.createDimension('lat',jpj)
        riv_a_n3n = ncfile.createVariable('riv_N3n', 'f4', ('lat','lon'))
        riv_a_n1p = ncfile.createVariable('riv_N1p', 'f4', ('lat','lon'))
        riv_a_n5s = ncfile.createVariable('riv_N5s', 'f4', ('lat','lon'))
        riv_a_o3c = ncfile.createVariable('riv_O3c', 'f4', ('lat','lon'))
        riv_a_o3h = ncfile.createVariable('riv_O3h', 'f4', ('lat','lon'))
        riv_a_O2o = ncfile.createVariable('riv_O2o', 'f4', ('lat','lon'))
        riv_a_Hg2 = ncfile.createVariable('riv_Hg2', 'f4', ('lat','lon'))
        riv_a_MMHg = ncfile.createVariable('riv_MMHg', 'f4', ('lat','lon'))
        riv_a_R3c = ncfile.createVariable('riv_R3c', 'f4', ('lat','lon'))
        riv_a_R3l = ncfile.createVariable('riv_R3l', 'f4', ('lat','lon'))

        riv_a_n3n[:] = self.get_map_from_1d_array(N, mask)
        riv_a_n1p[:] = self.get_map_from_1d_array(P, mask)
        riv_a_n5s[:] = self.get_map_from_1d_array(S, mask)
        riv_a_o3c[:] = self.get_map_from_1d_array(D, mask)
        riv_a_o3h[:] = self.get_map_from_1d_array(A, mask)
        riv_a_O2o[:] = self.get_map_from_1d_array(O, mask)
        riv_a_Hg2[:] = self.get_map_from_1d_array(H, mask)
        riv_a_MMHg[:] = self.get_map_from_1d_array(M, mask)
        riv_a_R3c[:] = self.get_map_from_1d_array(DOC, mask)
        riv_a_R3l[:] = self.get_map_from_1d_array(CDOM,mask)

        setattr(riv_a_n3n,'missing_value',np.float32(1.e+20))
        setattr(riv_a_n1p,'missing_value',np.float32(1.e+20))
        setattr(riv_a_n5s,'missing_value',np.float32(1.e+20))
        setattr(riv_a_o3c,'missing_value',np.float32(1.e+20))
        setattr(riv_a_o3h,'missing_value',np.float32(1.e+20))
        setattr(riv_a_O2o,'missing_value',np.float32(1.e+20))
        setattr(riv_a_Hg2,'missing_value',np.float32(1.e+20))
        setattr(riv_a_MMHg,'missing_value',np.float32(1.e+20))
        setattr(riv_a_R3c,'missing_value',np.float32(1.e+20))
        setattr(riv_a_R3l,'missing_value',np.float32(1.e+20))
        ncfile.close()
        return


    def get_map_from_1d_array(self,array,mask):
        '''
        Arguments:
        * array * an 1D numpuy array
        * mask  * Mask object
        We sum because there are more rivers that can be assigned to the same cell

        Returns:
        * OUT * a 2D map
        '''
        _, jpj, jpi = mask.shape
        OUT = np.ones((jpj, jpi), np.float32) * 1.e+20
        OUT[mask.mask_at_level(0)] = -1.
        if (self.nspread==1):
            for jr in range(self.nrivers):
                ji = self.georef['indLon'][jr] - 1
                jj = self.georef['indLat'][jr] - 1
                if OUT[jj,ji] < 0.:
                    OUT[jj,ji] = array[jr] # First time
                else:
                    OUT[jj,ji] += array[jr] # To handle cells with multiple river points
            return OUT
        else:
            for jr in range(self.nrivers*self.nspread):
                ji = self.georef_spread['indLon'][jr] - 1
                jj = self.georef_spread['indLat'][jr] - 1
                if OUT[jj,ji] < 0.:
                    OUT[jj,ji] = array[jr] # First time
                else:
                    OUT[jj,ji] += array[jr] # To handle cells with multiple river points
            return OUT


    def load_from_file(self,filename,var):
        ''' Useful for check/debug '''
        dset = netCDF4.Dataset(filename, 'r')
        M = np.array(dset[var])
        dset.close()
        return M
if __name__=="__main__":
    from bitsea.commons.mask import Mask
    import config as conf
    TheMask = Mask.from_file(conf.file_mask)
    R = river(conf)
    R.modularize(conf)
    R.gen_map_indexes(TheMask)
    
