#!/usr/bin/python3.4

import logging
from bclib.io_lib  import read_configure as rconf
from bclib.obj_lib import mask_obj as maskobj
from bclib.obj_lib import co2_obj as co2obj
from bclib.io_lib import excel_obj as xobj

logging.basicConfig( level=logging.DEBUG)
logging.info("BOUNDARY CONDITIONS python version 0.1")
logging.info(" 18 sett 2015 ")
logging.info(" logging set ")
elab = rconf.elaboration()

mesh = maskobj.mesh(elab)
mesh.generate_boundmask()
mesh.bounmesh.write_netcdf()
mesh.river.map_contribute_on_sea()
mesh.bc()
#mesh.submesh.atmosphere()
# co2  = co2obj.co2atm(elab.file_co2)
#mesh.generate_boundmask()
#mesh.bounmesh.write_netcdf()
# river = xobj.xlsx(elab.file_river)
# roundoff = xobj.xlsx(elab.file_runoff)

logging.info("end elaboration")