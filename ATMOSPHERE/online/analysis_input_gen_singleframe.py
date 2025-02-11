import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates daily files for OASIM by reading ECMWF file provided by CMCC
    Valid for analysis file having a single frame
    YYYYMMDD-ECMWF---AM0100-MEDATL-bYYYYMMDD_an00-fv12.00.nc analisi at 00 of the same day YYYYMMDD
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = ''' Input CMCC file'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' mask filename'''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' path of the output optical dir '''
                                )
    return parser.parse_args()


args = argument()


import os
from datetime import datetime, timedelta
from bitsea.commons import netcdf4
import numpy as np
from scipy import interpolate
from bitsea.commons.mask import Mask
from bitsea.commons.Timelist import TimeList
import netCDF4
from bitsea.commons.utils import addsep



TheMask=Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape

INPUTDIR=addsep(args.inputdir)
OUTDIR = addsep(args.outdir)

TL = TimeList.fromfilenames(None, INPUTDIR, "*_an00-fv12.00.nc", prefix="", dateformat="%Y%m%d")

nframes_in_day = 1
deltaH = 6

filename=TL.filelist[0]


lon= netcdf4.readfile(filename, 'lon')
lat= netcdf4.readfile(filename, 'lat')



xMin = TheMask.xlevels[0,0]
xMax = TheMask.xlevels[0,-1]
yMin = TheMask.ylevels[0,0]
yMax = TheMask.ylevels[-1,0]


I_start = np.argmin(np.abs(lon - xMin))
I_end   = np.argmin(np.abs(lon - xMax))
J_start = np.argmin(np.abs(lat - yMax))
J_end   = np.argmin(np.abs(lat - yMin))

lon = lon[I_start:I_end]
lat = lat[J_start:J_end]
lat = lat[-1::-1] # reverse


def readframe(filename,var, timeframe):
    A=netcdf4.readfile(filename, var)[timeframe,J_start:J_end,I_start:I_end]
    return A[-1::-1,:]

def interp(Min):
    X,Y = np.meshgrid(lon,lat)
    nPoints = X.size
    Xpoints = np.zeros((nPoints,2),float)
    Xpoints[:,0] = X.ravel()
    Xpoints[:,1] = Y.ravel()
    f = interpolate.NearestNDInterpolator(Xpoints, Min.ravel())
    Mout  = f(TheMask.xlevels, TheMask.ylevels)
    return Mout

def getframe(filename,var, timeframe):
    M2d_orig = readframe(filename, var, timeframe)
    return  interp(M2d_orig)

def dumpfile(filename, maskObj, sp,msl, t2m,d2m, tcc,w10):
    ncOUT   = netCDF4.Dataset(filename,"w");

    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);

    ncvar=ncOUT.createVariable('lon','f',('lon',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'longitude')
    ncvar[:]=maskObj.xlevels[0,:]
    ncvar=ncOUT.createVariable('lat','f',('lat',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'latitude')
    ncvar[:]=maskObj.ylevels[:,0]


    ncvar = ncOUT.createVariable('sp','f',('lat','lon'))
    ncvar[:]=sp
    setattr(ncvar, 'long_name',  'Surface pressure' )
    setattr(ncvar, 'units','Pa' )
    setattr(ncvar, 'code', 34 )
    

    ncvar = ncOUT.createVariable('msl','f',('lat','lon'))
    setattr(ncvar, 'units','Pa' )
    setattr(ncvar,'long_name','Mean sea-level pressure')
    setattr(ncvar,'code',151)
    ncvar[:] = msl

    ncvar = ncOUT.createVariable('d2m','f',('lat','lon'))
    setattr(ncvar,'units','K')
    setattr(ncvar,'long_name','2 metre dewpoint temperature')
    setattr(ncvar, 'code', 168)
    ncvar[:] = d2m

    ncvar = ncOUT.createVariable('t2m','f',('lat','lon'))
    setattr(ncvar,'units','K')
    setattr(ncvar,'long_name','2 metre temperature')
    setattr(ncvar, 'code', 167)
    ncvar[:] = t2m

    ncvar = ncOUT.createVariable('tcc','f',('lat','lon'))
    setattr(ncvar,'units','%')
    setattr(ncvar,'long_name','Total cloud cover')
    setattr(ncvar, 'code', 164)
    ncvar[:] = tcc

    ncvar = ncOUT.createVariable('w10','f',('lat','lon'))
    setattr(ncvar,'units','m/s')
    setattr(ncvar,'long_name','10 metre wind speed module')
    setattr(ncvar,'code', '165 and 166')
    ncvar[:] = w10

    setattr(ncOUT, 'input_file', str(inputfile))
    ncOUT.close()


for inputfile in TL.filelist:
    yyyymmdd=os.path.basename(inputfile)[:8]
    for iframe in range(nframes_in_day):
        d=datetime.strptime(yyyymmdd,'%Y%m%d') + timedelta(hours = (deltaH*iframe + deltaH/2))
        outfile = OUTDIR + d.strftime("atm.%Y%m%d-%H:%M:%S.nc")
        print(outfile)


        msl = getframe(inputfile,'MSL' , iframe)
        sp =  getframe(inputfile,'SP'  , iframe)
        u10 = getframe(inputfile,'U10M', iframe)
        v10 = getframe(inputfile,'V10M', iframe)
        t2m = getframe(inputfile,'T2M' , iframe)
        d2m = getframe(inputfile,'D2M' , iframe)
        tcc = getframe(inputfile,'TCC' , iframe)
    
        w10 = np.sqrt(u10**2 + v10**2)
        dumpfile(outfile, TheMask, sp,msl, t2m,d2m, tcc,w10)
    
