import numpy as np
from bitsea.commons import netcdf4
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
import seawater as sw
import netCDF4
from bclib.river import conversion


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


INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/PO/TEST_01/wrkdir/MODEL/FORCINGS/"
OUTDIR="/gpfs/scratch/userexternal/gbolzon0/PO/ogstm_boundary_conditions/out/"

riverinput="/gpfs/scratch/userexternal/gbolzon0/PO/runoff_1d_nomask_y2019.nc"

CMCC_Mask=Mask.from_file('/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask_CMCCfor_ogstm.nc')
TheMask=Mask.from_file('/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask.nc')
jpk, jpj, JPI = CMCC_Mask.shape
jpk, jpj, jpi = TheMask.shape
mask0 = TheMask.mask_at_level(0)

S = netcdf4.readfile(riverinput, 's_river')

ii = S[0,:,:]==18
J,I = np.nonzero(ii)

VARLIST=['O2o','N1p','N3n','N5s','O3c','O3h']

RIVER_CONCENTRATION={'O2o':250  ,  #mmol/m3
                     'N1p':2.572,  #mmol/m3
                     'N3n':150  ,  #mmol/m3
                     'N5s':75   ,  #mmol/m3
                     'O3c':33225,  #  mg/m3
                     'O3h':2800 }  #mmol/m3



nPoints=ii.sum()
SALI = np.ones((10,),np.float32)*18
PRES = np.ones((10,),np.float32)*CMCC_Mask.zlevels[0]

cut =JPI - jpi
TI = TimeInterval('2019','2020','%Y')
TL = TimeList.fromfilenames(TI, INPUTDIR, "T*.nc", prefix="T")
nFrames = TL.nTimes
my_dtype=[(var, np.float32) for var in VARLIST]
M2D = np.ones((jpj,JPI), dtype=my_dtype)
PO  = np.ones((jpj,jpi), dtype=my_dtype)
Total_Po = np.zeros((nFrames,), dtype=my_dtype)

yearly_volume=0

for iFrame, filename in enumerate(TL.filelist):
    outfile=OUTDIR + "PO__" + TL.Timelist[iFrame].strftime("%Y%m%d-%H:%M:%S") + ".nc"
    print outfile
    Runoff   = DataExtractor(CMCC_Mask, filename, 'sorunoff', dimvar=2).values
    TEMP     = DataExtractor(CMCC_Mask, filename, 'votemper').values[0,:]

    T      = sw.temp(SALI,TEMP[ii],PRES)
    RHO    = sw.dens(SALI,T,PRES)
    
    Discharge = Runoff[ii]*CMCC_Mask.area[ii]/RHO  # m^3/s
    daily_volume = Discharge.sum()*86400
    yearly_volume += daily_volume

    
    for var in VARLIST:
        Total_Po[iFrame][var] = daily_volume * RIVER_CONCENTRATION[var] # mmol/day
        M2D[var] = -1.0
        M2D[var][ii] = Discharge * RIVER_CONCENTRATION[var] / CMCC_Mask.area[ii]
        
        PO[var] = M2D[var][:,cut:]
        PO[var][~mask0] = 1.e+20
        
    dump_file(outfile, PO['N3n'], PO['N1p'], PO['N5s'], PO['O3h'], PO['O3c'], PO['O2o'], TheMask)


np.save('Total_po.npy',Total_Po)


for var in VARLIST:
    print "mmol/year", var, yearly_volume*RIVER_CONCENTRATION[var]
    print "KTons/year", var, conversion(var)*yearly_volume*RIVER_CONCENTRATION[var]


# KTons/year O2o 12.215648423583986
# KTons/year N1p 3.8959123204367927
# KTons/year N3n 102.61144675810547
# KTons/year N5s 205.22289351621095
# KTons/year O3c 1623.4596754943116
# KTons/year O3h 136.81526234414062
