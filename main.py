#!/usr/bin/python3.5
from bclib.io_lib  import read_configure as rconf
from bclib.obj_lib import mask_obj as maskobj


elab = rconf.elaboration(json_input="./conf24.json")


mesh = maskobj.mesh(elab)
#mesh.generate_boundmask()
#mesh.bounmesh.write_netcdf()
#mesh.submesh.atmosphere()
#mesh.river.map_contribute_on_sea()
#mesh.bc()
#co2  = co2obj.co2atm(elab)
#co2.generator(mesh)
