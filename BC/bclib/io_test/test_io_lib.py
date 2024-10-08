import numpy as np
from bclib.io_lib import read_configure as rconf

def test_json():
    """ Test json reading class, the input is a string in this case """

    example = """{
        "nameElaboration": "prova1",
        "parameter" :
         {
          "simulationTime":1000,
          "start_time":1960,
          "end_time":2000,
          "num_variables" : 7,
          "CO2_end" : 337,
          "CO2_start" : 136
         },

        "filepath" :
        {
        "file_river":"/home/eric/ogs_tools/matfunc/bc/INPUT_DATA/TI/D4.3.2_SESAME/N_P_SES_rivers.xlsx",
        "file_runoff":"/home/eric/ogs_tools/matfunc/bc//INPUT_DATA/TI/D4.3.2_SESAME/N_P_SES_runoff.xlsx",
        "file_nutrients":"/home/eric/ogs_tools/matfunc/bc//INPUT_DATA/GIB/VPnutrients_CO2.nc",
        "file_co2":"/home/eric/ogs_tools/matfunc/bc/INPUT_DATA/CO2/CMIP5_scenarios_RCP_CO2_mixing_ratio.nc",
        "dir_out":"./"
        },

        "atmosphere" :
        {
          "N3n_wes":41433.5,
          "N3n_eas":39842.5,
          "P04_wes":527,
          "P04_eas":668
        },
        
        "river" :
        {
          "sheet_array":
          [
            {"name":"km3_per_yr"},
            {"name":"no3_kt_yr"},
            {"name":"po4_kt_yr"},
            {"name":"dic_kt_yr"},
            {"name":"alk_Gmol_yr"}

          ]
        },

        "variables" :
        {
          "end_nudging":-5.25,
          "rdpmax":0.5,
          "rdpmin":0.0416,
          "var_array":
          [
            {"name":"N1p","end_max":-7.5},
            {"name":"N3n","end_max":-7.5},
            {"name":"O2o","end_max":-7.5},
            {"name":"N5s","end_max":-7.5},
            {"name":"O3h","end_max":-5.5},
            {"name":"O3h","end_max":-5.5},
            {"name":"N6r","end_max":-5.5}
          ]
       },



        "mask":
          {
            "nameMask" : "v1",
            "maskfile" :  "/home/eric/ogs_tools_input/meshmask_872.nc",
            "submaskfile"  : "/home/eric/ogs_tools_input/sbmask_872.nc",
            "bounmask" : "/home/eric/ogs_tools_input/sbmask_872.nc",
            "jpi" : 362,
            "jpj" : 128,
            "jpk" : 72,
            "lon_start" :  -8.78125,
            "lon_step" :  0.125,
            "lat_start" :  30.15620,
            "lat_step" :  0.125
           }


    }"""

    json_rid = rconf.elaboration(example)
    assert json_rid.simulation_start_time == 1960
    assert json_rid.variables[0][0] == "N1p"
    assert json_rid.variables[0][1] == -7.5

    test_lon_np_array = np.arange(-8.78125,-8.78125+362*0.125,0.125)
    test_lat_np_array = np.arange(30.15620,30.15620+128*0.125,0.125)
    datasheet = ["km3_per_yr","no3_kt_yr","po4_kt_yr","dic_kt_yr","alk_Gmol_yr"]
    
    assert (json_rid.lon_np_array == test_lon_np_array).all()
    assert (json_rid.lat_np_array == test_lat_np_array).all()
    assert json_rid.river_data_sheet == datasheet
