import numpy as np
import scipy.io.netcdf as nc
from bclib.obj_lib import mask_obj
import datetime
from matplotlib.dates import seconds
from scipy.constants.constants import atmosphere
from bclib.obj_lib.mask_obj import sub_mesh
#from ctypes.wintypes import BOOLEAN

class co2atm:
    """
        Class co2

    """

    def __init__(self, input_data):
        self.input_data = input_data
        self.path = input_data.file_co2
        self._extract_information()


    def _extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()
    
    def generator(self, mesh):
        time = self.time[self.input_data.co2_start:self.input_data.co2_end]
        rcp45 = self.RCP45[self.input_data.co2_start:self.input_data.co2_end]
        rcp85 = self.RCP85[self.input_data.co2_start:self.input_data.co2_end]
        
        date_format = "%Y %m %d - %H:%M:%S"                
               
        factor = 60*60*24  #factor of conversion seconds to days
        giorni = []
        nb = 0             #number of bisestile years from 0 to 1765
        nnb = 0            #number of not bisestile years from 0 to 1765
        
        for i in range(0, 1766):        #calculation of bisestile and non bisestile years from 0 to 1765
            if divmod(i,4)[1] ==0:
                nb +=1
            else:
                nnb+=1
        
        
        
        total_days = 365*nnb + 366*nb  #total days from year 0 to 1765
        
        for sec in time:    
            giorni.append(total_days + sec/factor)
        
        #print(total_days, nnb, nb) it was just a check of the numbers
                    
        co2datenumber = []
        for d in giorni:
            #if long(d) <= datetime.datetime.max.toordinal():
            co2datenumber.append(datetime.datetime.fromordinal(int(d)))
              
        co2datestr = []
        for dn in co2datenumber:
            co2datestr.append(dn.strftime(date_format))
        
        
        
        count = 0
        
        #l_tmask = bool(mesh.tmask.all()) WRONG!!!
        l_tmask = mesh.tmask.astype(np.bool)


        for yCO2 in co2datestr:
        
            fileOUT = self.input_data.dir_out + "/CO2_" + yCO2 + ".nc"
            map_co2 = np.dot(np.ones([self.input_data.jpj, self.input_data.jpi]), rcp85[count])
            for ciccio in l_tmask:
                map_co2[not(ciccio)] = np.nan
            ncfile = nc.netcdf_file(fileOUT, 'w')
            ncfile.createDimension('lon', self.input_data.jpi)
            ncfile.createDimension('lat', self.input_data.jpj)
            g = ncfile.createVariable('co2', 'f', ('lon','lat'))
            g = map_co2
            #g.longname = "CO2 content" WRONG!
            setattr(self, 'longname', "CO2 content")
            setattr(ncfile, 'date', datetime.datetime.now().strftime(date_format))
            setattr(ncfile, 'comment', "Uniform value")
            ncfile.close()
            count +=1
            
            
        