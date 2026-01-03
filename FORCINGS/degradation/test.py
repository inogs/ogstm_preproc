import xarray as xr
from degrade_mesh import xpnd_wrap, degr_wrap
from degrade_mesh import jsum_istep
M = xr.open_dataset(maskfile)
ndeg=6
a = xpnd_wrap(M["e2t"], 'edge', ndeg) # in load_mesh
e2t = degr_wrap(a, jsum_istep, ndeg, W=None)
