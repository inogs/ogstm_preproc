from config import ICdef, LayerList
from commons.mask import Mask
from climatology import get_climatology
import numpy as np
import pylab as pl
from commons.submask import SubMask
import scipy.io.netcdf as NC
import density


maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc"
maskfile_ingv="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask_sizeINGV.nc"
filename="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/FORCINGS/eas2_v12_notint_1d_20140101_20140102_grid_T.nc"

TheMask = Mask(maskfile)
IngvMask= Mask(maskfile_ingv)

PresCEN = np.array([(l.bottom+l.top)/2  for l in LayerList])
jpk, jpj, jpi = TheMask.shape
JPK, JPJ, JPI = IngvMask.shape

def shift(M2d,pos, verso):
    out = np.ones_like(M2d)*np.nan
    
    if (verso=='u'): out[:-pos,:] = M2d [pos:,:]
    if (verso=='d'): out[ pos:,:] = M2d[:-pos,:]
    if (verso=='l'): out[:,:-pos] = M2d [:,pos:]
    if (verso=='r'): out[:,pos:]  = M2d[:,:-pos]

    if (verso=="ul") : out[:-pos,:-pos] =  M2d [pos:,pos:]
    if (verso=="ur") : out[:-pos,pos: ] =  M2d [pos:,:-pos]
    if (verso=="dl") : out[ pos:,pos:]  =  M2d [:-pos,pos:]
    if (verso=="dr") : out[pos:,pos: ]  = M2d [:-pos,:-pos]
    return out
    


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
    ncOUT=NC.netcdf_file(outfile,"w")
    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',1)

    TRN   = 'TRN' + var;
    ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = TheMask.ylevels
    ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = TheMask.zlevels
    ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
    ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x')); ncvar[:] = rst;


    setattr(ncOUT.variables[TRN]   ,'missing_value',1e+20                              );
    setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
    setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
    ncOUT.close()

    

rho =density.get_density(filename, IngvMask)
rho=rho[:jpk,:,JPI-jpi:]

VARLIST=['N1p','N3n','O2o','N5s','O3h','O3c']

for varname in VARLIST:
    CLIM = get_climatology(varname)
    out_img = varname + ".png"
    outfile = "RST.20140101-00:00:00." + varname + ".nc"
    print outfile
    RST = np.zeros((jpk,jpj,jpi),np.double)
     
    fig, ax = pl.subplots()
    for isub, sub in enumerate(ICdef):
        p = getModelProfile(CLIM[isub,:])
        ax.plot(p, TheMask.zlevels,'b')
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
    ax.invert_yaxis()
    fig.savefig(out_img)


