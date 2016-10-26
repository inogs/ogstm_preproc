#!/usr/bin/python3.5
from bclib.io_lib  import read_configure as rconf
#from bclib.obj_lib import mask_obj as maskobj
#from commons.mask import Mask

conf = rconf.elaboration(json_input="./conf24.json")

conf.file_river= "input_obc_eas2_v3.xlsx"
from bclib.obj_lib.river_obj import river_data
R = river_data(conf)
R.modularize(conf)
R.gen_map_indexes(conf)
import sys
sys.exit()

# for sheet in conf.river_data_sheet:
#     R.river_data[sheet]['yyyy'] = (R.river_data[sheet]['2000'] + R.river_data[sheet]['2001'] + R.river_data[sheet]['2002'])/3
# R.generate_climatological_monthly_files(conf, mask)


#mask = Mask(conf.file_mask)

#mesh = maskobj.mesh(elab)
#mesh.generate_boundmask()
#mesh.bounmesh.write_netcdf()
#mesh.submesh.atmosphere()
#mesh.river.map_contribute_on_sea()
#mesh.bc()
#co2  = co2obj.co2atm(elab)
#co2.generator(mesh)
