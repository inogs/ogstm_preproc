from config import ICdef, LayerList
from commons.mask import Mask
from climatology import get_climatology
import numpy as np
import pylab as pl
from commons.submask import SubMask
from IC import RSTwriter, smoother
from commons import density
from commons.utils import getcolor



def getModelProfile(climvalues):
    zclim = PresCEN
    ii=np.isnan(climvalues)
    zclim = zclim[~ii]
    climvalues = climvalues[~ii]
    return np.interp(TheMask.zlevels,zclim,climvalues)





maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc"
maskfile_ingv="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask_sizeINGV.nc"
filename="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/FORCINGS/eas2_v12_notint_1d_20140101_20140102_grid_T.nc"

TheMask = Mask(maskfile)
IngvMask= Mask(maskfile_ingv)

PresCEN = np.array([(l.bottom+l.top)/2  for l in LayerList])
jpk, jpj, jpi = TheMask.shape
JPK, JPJ, JPI = IngvMask.shape
rho =density.get_density(filename, IngvMask)
rho=rho[:jpk,:,JPI-jpi:]

VARLIST=['N1p','N3n','O2o','N5s','O3h','O3c']


nSub = len(ICdef.basin_list)
# first of all, plots
for varname in VARLIST:
    CLIM = get_climatology(varname)
    out_img = varname + ".png"
    fig, ax = pl.subplots()
    for isub, sub in enumerate(ICdef):
        p = getModelProfile(CLIM[isub,:])
        ax.plot(p, TheMask.zlevels,color=getcolor(nSub, isub), label=sub.name)
    ax.invert_yaxis()
    ax.legend()
    fig.savefig(out_img)
    pl.close(fig)

# then 3D arrays generation ---------
for varname in VARLIST:
    CLIM = get_climatology(varname)
    outfile = "RST.20140101-00:00:00." + varname + ".nc"
    print outfile
    RST = np.zeros((jpk,jpj,jpi),np.double)
    for isub, sub in enumerate(ICdef):
        p = getModelProfile(CLIM[isub,:])
        S =SubMask(sub, maskobject=TheMask)
        for k in range(jpk):
            submask = S.mask[k,:,:]
            V = RST[k,:,:]
            V[submask] =p[k]
    if varname == 'O3h':  RST = RST*rho/1000
    if varname == 'O3c':  RST = RST*rho*12/1000
    print "smoother"
    RST_s = smoother(TheMask, RST)
    print "writer"
    check = np.isnan(RST_s[TheMask.mask])
    print "number of nans: ", check.sum()
    RSTwriter(outfile, varname, RST_s, TheMask)

