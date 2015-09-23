import json
import numpy as np
import os.path


class elaboration:

    """
    This function read the configuration file that set the
    elaboration parameter.

    json_input(optional) -> path of file that contain configuration,
                          the default path is ./conf.json

    """

    def __init__(self, json_input="./conf.json"):
        """ Constructor, read json string and call set function """

        if os.path.isfile(json_input) :
            self.init_configure = json.loads(open(json_input).read())
        else :
            self.init_configure = json.loads(json_input)
            print("stringa")


        self.read_parameter_section()
        self.read_filepath_section()
        self.read_atmosphere()
        self.read_variable()
        self.read_mask()

    def read_parameter_section(self):
        """ Set general elaboration parameter """
        self.name = self.init_configure["nameElaboration"]
        param = self.init_configure["parameter"]
        # for i in param:
        #     setattr(self, i, param[i])
        self.simulationTime = param["simulationTime"]
        self.simulation_start_time = param["start_time"]
        self.simulation_end_time = param["end_time"]
        self.co2_start = param["CO2_start"]
        self.co2_end = param["CO2_end"]

    def read_filepath_section(self):
        """ Set general file path parameter """
        file_path = self.init_configure["filepath"]
        # for i in file_path:
        #     setattr(self, i, file_path[i])
        self.file_river = file_path["file_river"]
        self.file_runoff = file_path["file_runoff"]
        self.file_nutrients = file_path["file_nutrients"]
        self.file_co2 = file_path["file_co2"]
        self.dir_out = file_path["dir_out"]

    def read_atmosphere(self):
        """ Set general atmosphere parameter """
        atm = self.init_configure["atmosphere"]
        self.n3n_wes = atm["N3n_wes"]
        self.n3n_eas = atm["N3n_eas"]
        self.po4_wes = atm["P04_wes"]
        self.po4_eas = atm["P04_eas"]

    def read_variable(self):
        """ Set variables and their parameter """
        var = self.init_configure["variables"]
        self.rdpmax = var["rdpmax"]
        self.rdpmin = var["rdpmin"]
        self.end_nudging = var["end_nudging"]
        self.variables = [ ]
        for i in var["var_array"]:
            self.variables.append([i["name"],i["end_max"]])

    def read_mask(self):
        """ Set mask parameter """
        msk =  self.init_configure["mask"]
        self.name_mask = msk["nameMask"]
        self.file_mask = msk["maskfile"]
        self.file_submask = msk["submaskfile"]
        self.file_bmask = msk["bounmask"]
        self.jpi = msk["jpi"]
        self.jpj = msk["jpj"]
        self.jpk = msk["jpk"]
        self.lon_np_array = np.arange(msk["lon_start"],msk["lon_start"]+msk["jpi"]*msk["lon_step"],msk["lon_step"])
        self.lat_np_array = np.arange(msk["lat_start"],msk["lat_start"]+msk["jpj"]*msk["lat_step"],msk["lat_step"])
