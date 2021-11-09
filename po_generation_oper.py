import argparse


def argument():
    parser = argparse.ArgumentParser(description = '''Generates PO files reading runoff from forcings
                               ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help ='/some/path/'
                                )
    parser.add_argument(   '--outdir',"-o",
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--cmccmaskfile', '-M',
                                type = str,
                                required = True,
                                help = '/some/path/outmask.nc')
    parser.add_argument(   '--maskfile',
                                type = str,
                                required = True,
                                help = '''NetCDF File name of the mask.
                                ''')
    return parser.parse_args()


args = argument()
from commons.utils import addsep
import numpy as np
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.dataextractor import DataExtractor
import seawater as sw
import netCDF4


def dump_file(filename,N,P,S,A,D,O,mask):
    '''
      Writes the single PO_ file
    '''
    _,jpj, jpi = mask.shape
    ncfile = netCDF4.Dataset(filename, 'w')
    ncfile.createDimension('lon',jpi)
    ncfile.createDimension('lat',jpj)
    riv_a_n3n = ncfile.createVariable('riv_N3n', 'f4', ('lat','lon'))
    riv_a_n1p = ncfile.createVariable('riv_N1p', 'f4', ('lat','lon'))
    riv_a_n5s = ncfile.createVariable('riv_N5s', 'f4', ('lat','lon'))
    riv_a_o3c = ncfile.createVariable('riv_O3c', 'f4', ('lat','lon'))
    riv_a_o3h = ncfile.createVariable('riv_O3h', 'f4', ('lat','lon'))
    riv_a_O2o = ncfile.createVariable('riv_O2o', 'f4', ('lat','lon'))
    riv_a_n3n[:] = N
    riv_a_n1p[:] = P
    riv_a_n5s[:] = S
    riv_a_o3c[:] = D
    riv_a_o3h[:] = A
    riv_a_O2o[:] = O

    setattr(riv_a_n3n,'missing_value',np.float32(1.e+20))
    setattr(riv_a_n1p,'missing_value',np.float32(1.e+20))
    setattr(riv_a_n5s,'missing_value',np.float32(1.e+20))
    setattr(riv_a_o3c,'missing_value',np.float32(1.e+20))
    setattr(riv_a_o3h,'missing_value',np.float32(1.e+20))
    setattr(riv_a_O2o,'missing_value',np.float32(1.e+20))
    ncfile.close()
    return

INPUTDIR=addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
CMCC_Mask=Mask(args.cmccmaskfile)
TheMask=Mask(args.maskfile)
jpk, jpj, JPI = CMCC_Mask.shape
jpk, jpj, jpi = TheMask.shape
mask0 = TheMask.mask_at_level(0)


J=np.array([350, 350, 350, 352, 353, 354, 355, 356, 357, 358])
I=np.array([730, 732, 733, 735, 736, 737, 737, 735, 733, 732])
VARLIST=['O2o','N1p','N3n','N5s','O3c','O3h']

RIVER_CONCENTRATION={'O2o':250  ,  #mmol/m3
                     'N1p':2.572,  #mmol/m3
                     'N3n':150  ,  #mmol/m3
                     'N5s':75   ,  #mmol/m3
                     'O3c':33225,  #  mg/m3
                     'O3h':2800 }  #mmol/m3



nPoints=len(I)
SALI = np.ones((nPoints,),np.float32)*18
PRES = np.ones((nPoints,),np.float32)*CMCC_Mask.zlevels[0]
cut =JPI - jpi

my_dtype=[(var, np.float32) for var in VARLIST]
M2D = np.ones((jpj,JPI), dtype=my_dtype)
PO  = np.ones((jpj,jpi), dtype=my_dtype)
TL = TimeList.fromfilenames(None, INPUTDIR, "T*.nc", prefix="T")

ii = np.zeros((jpj,JPI),np.bool)
for k in range(nPoints): 
    ii[J[k],I[k]]=True


for iFrame, filename in enumerate(TL.filelist):
    outfile=OUTDIR + "PO__" + TL.Timelist[iFrame].strftime("%Y%m%d-%H:%M:%S") + ".nc"
    print outfile
    Runoff   = DataExtractor(CMCC_Mask, filename, 'sorunoff', dimvar=2).values
    TEMP     = DataExtractor(CMCC_Mask, filename, 'votemper').values[0,:]

    T      = sw.temp(SALI,TEMP[ii],PRES)
    RHO    = sw.dens(SALI,T,PRES)
    
    Discharge = Runoff[ii]*CMCC_Mask.area[ii]/RHO  # m^3/s
    #daily_volume = Discharge.sum()*86400

    
    for var in VARLIST:
        M2D[var] = -1.0
        M2D[var][ii] = Discharge * RIVER_CONCENTRATION[var] / CMCC_Mask.area[ii]
        
        PO[var] = M2D[var][:,cut:]
        PO[var][~mask0] = 1.e+20
        
    dump_file(outfile, PO['N3n'], PO['N1p'], PO['N5s'], PO['O3h'], PO['O3c'], PO['O2o'], TheMask)






