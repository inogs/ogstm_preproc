#!/usr/bin/python3.5
from bclib.io_lib  import read_configure as rconf
from bclib.obj_lib.bmask_obj import bmesh
#from bclib.obj_lib import mask_obj as maskobj
from commons.mask import Mask

conf = rconf.elaboration(json_input="./conf24.json")

TheMask = Mask(conf.file_mask)
B=bmesh(conf.file_bmask, conf)
index = B.load_bounmask()

conf.file_river= "input_obc_eas2_v3.xlsx"
from bclib.obj_lib.river_obj import river_data
R = river_data(conf)
R.modularize(conf)
georef = R.gen_map_indexes(TheMask)
idxt, positions = R.gen_boun_indexes(index)
import sys
sys.exit()



R.generate_monthly_files(conf, TheMask,idxt, positions)

climatological=False
if climatological:
    for sheet in conf.river_data_sheet:
        R.river_data[sheet]['yyyy'] = (R.river_data[sheet]['2000'] + R.river_data[sheet]['2001'] + R.river_data[sheet]['2002'])/3
    R.generate_climatological_monthly_files(conf, TheMask)




#mesh = maskobj.mesh(elab)
#mesh.generate_boundmask()
#mesh.bounmesh.write_netcdf()
#mesh.submesh.atmosphere()
#mesh.river.map_contribute_on_sea()
#mesh.bc()
#co2  = co2obj.co2atm(elab)
#co2.generator(mesh)
