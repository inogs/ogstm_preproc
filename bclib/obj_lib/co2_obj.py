
import numpy as np
import scipy.io.netcdf as nc

class co2atm:
    """
        Class co2

    """

    def __init__(self, ncfile):
        self.path = ncfile
        self.extract_information()


    def extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
             setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            setattr(self, i, self.ncfile.variables[i])
