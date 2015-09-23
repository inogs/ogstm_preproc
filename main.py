#!/usr/bin/python3.4


from bclib.io_lib import read_configure as rconf
from bclib.obj_lib import mask_obj as maskobj
from bclib.obj_lib import co2_obj as co2obj
from bclib.obj_lib import excel_obj as xobj


print("BOUNDARY CONDITIONS python version 0.1")
print(" 18 sett 2015 ")

elab = rconf.elaboration()
mesh = maskobj.mesh(elab.file_mask,elab.file_submask,elab.file_bmask)
co2  = co2obj.co2atm(elab.file_co2)
mesh.generate_boundmask(elab)
mesh.bounmesh.write_netcdf()
# river = xobj.xlsx(elab.file_river)
# roundoff = xobj.xlsx(elab.file_runoff)

print("FINISH")