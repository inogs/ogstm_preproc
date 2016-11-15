import numpy as np
import netCDF4 as nc
import datetime
from commons.time_interval import TimeInterval

import logging

class co2atm():
    """
        Class co2atm
        author : Giorgio Bolzon
    """

    def __init__(self, conf):
        logging.info("CO2 enabled")
        self.config = conf
        self.path = conf.file_co2
        self._extract_information()
        basetime = datetime.datetime(1765,1,1,0,0,0)
        timelist=[]
        for sec in self.time:
            timeobj = basetime + datetime.timedelta(seconds= sec)
            timelist.append(timeobj)
        self.timelist = timelist


    def _extract_information(self):
        try:
            self.ncfile = nc.Dataset(self.path, 'r')
        except:
            print("CO2 FILE NOT FOUND")
            exit()

        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()

    def generate(self, mask):

        l_tmask = ~ mask.mask_at_level(0)
        _, jpj, jpi = mask.shape

        starttime = str(self.config.co2_start-1) + "0101"
        end__time = str(self.config.co2_end  +2) + "0101"
        TI = TimeInterval(starttime,end__time,"%Y%m%d")

        for it, t in enumerate(self.timelist):
            if TI.contains(t):
                fileOUT = self.config.dir_out + "/CO2_" + t.strftime("%Y%m%d") + "-12:00:00.nc"
                map_co2 = np.ones((jpj,jpi), dtype=np.float32) *self.RCP85[it]
                map_co2[l_tmask] = 1.e+20

                ncfile = nc.Dataset(fileOUT, 'w')
                ncfile.createDimension('lon', jpi)
                ncfile.createDimension('lat', jpj)
                ncvar = ncfile.createVariable('co2','f', ('lat','lon'))
                ncvar[:] = map_co2[:]
                setattr(ncvar,"missing_value",1e+20)
                setattr(self, 'longname', "CO2 content")
                setattr(ncfile, 'date', datetime.datetime.now().strftime("%Y%m%d-%H:%M:%S"))
                setattr(ncfile, 'comment', "Uniform value")
                ncfile.close()


        logging.info("CO2 files written")
if __name__ == "__main__":
    from bclib.io_lib  import read_configure
    from commons.mask import Mask
    conf = read_configure.elaboration(json_input="../../conf24.json")
    conf.file_mask="../../masks/meshmask_872.nc"
    conf.dir_out = "../../out"
    TheMask = Mask(conf.file_mask)
    conf.file_co2="../../CMIP5_scenarios_RCP_CO2_mixing_ratio.nc"
    CO2 =co2atm(conf)
    CO2.generate(TheMask)
