from config import ICdef, LayerList
from commons.mask import Mask
from climatology import get_climatology
import numpy as np
import pylab as pl
from commons.submask import SubMask
import netCDF4
import density
from commons.utils import getcolor
from commons.interpolators import shift


def getModelProfile(climvalues):
    zclim = PresCEN
    ii=np.isnan(climvalues)
    zclim = zclim[~ii]
    climvalues = climvalues[~ii]
    return np.interp(TheMask.zlevels,zclim,climvalues)

def smoother(maskobj,RST):
    NsmoothingCycles = 20
    jpk,jpj,jpi = maskobj.shape
    RST[~maskobj.mask] = np.nan
    for k in range(jpk):
        var2d = RST[k,:,:]
        auxs = np.zeros((5,jpj,jpi),np.float32)
        for _ in range(NsmoothingCycles):
            auxs[0,:,:] = var2d
            auxs[1,:,:] = shift(var2d,1,'r')
            auxs[2,:,:] = shift(var2d,1,'l')
            auxs[3,:,:] = shift(var2d,1,'u')
            auxs[4,:,:] = shift(var2d,1,'l')
            var2d = np.nanmean(auxs,axis=0)
        RST[k,:,:] = var2d
    return RST
    

def RSTwriter(outfile, var,rst):
    rst[~TheMask.mask] = 1.e+20
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',1)

    TRN   = 'TRN' + var;
    ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = TheMask.ylevels
    ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = TheMask.zlevels
    ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
    ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x'),zlib=True); ncvar[:] = rst;


    setattr(ncOUT.variables[TRN]   ,'missing_value',1e+20                              );
    setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
    setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
    ncOUT.close()


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
    print "smooter"
    RST_s = smoother(TheMask, RST)
    print "writer"
    check = np.isnan(RST_s[TheMask.mask])
    print "number of nans: ", check.sum()
    RSTwriter(outfile, varname, RST_s)

