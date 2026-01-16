# to load them it is needed the following virtual environment:
# source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh

import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates a restart for R1l, R2l R3l, R3c
    ''')
    parser.add_argument(   '--outfile_prefix', '-o',
                                type = str,
                                required = True,
                                default = "out/dir/RST.19950101-00:00:00",
                                help = 'Path of outfiles without .R1l.nc'
                                )
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )    
    return parser.parse_args()

args = argument()



from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.basins import V2
from bitsea.basins.basin import ComposedBasin
from IC import RSTwriter
import numpy as np
from bitsea.commons import netcdf4


TheMask=Mask.from_file(args.maskfile)

def CDOM_profile(surface_value,bottom_value,depth,depth0):
    k=0.05
    b=bottom_value/surface_value -1.0
    res=surface_value + (bottom_value-surface_value)/(1.0+np.exp(-k*(depth-depth0)))
    return res

EAS=ComposedBasin('eas',[V2.eas,V2.adr1,V2.adr2,V2.aeg],'East with marginal seas')
# Here we use basin objet features to run faster
EAS.region.border_longitudes    = [ 9,36]
EAS.region.border_latitudes     = [30,46]
V2.wes.region.border_longitudes = [-6,18]
V2.wes.region.border_latitudes  = [32,45]
#---------------------------------------------
print("wes")
wes = SubMask(V2.wes,mask=TheMask)
print("eas")
eas = SubMask(EAS,   mask=TheMask)


jpk,jpj,jpi=TheMask.shape
#R1l = np.ones((jpk,jpj,jpi), np.float64)*0.0
#R2l = np.ones((jpk,jpj,jpi), np.float64)*0.0
R3l = np.ones((jpk,jpj,jpi), np.float64)*1.0

nav_lev=TheMask.zlevels
for jk,depth in enumerate(nav_lev):
    print(depth)
    print('profile value west: ', CDOM_profile(0.8,1.0,depth,150.))
    print('profile value east: ', CDOM_profile(0.4,0.9,depth,150.))
    R3l[jk,wes.mask[jk,:,:]]=CDOM_profile(0.8,1.0,depth,150.)
    R3l[jk,eas.mask[jk,:,:]]=CDOM_profile(0.4,0.9,depth,150.)

#RSTwriter(args.outfile_prefix + ".R1l.nc", "R1l", R1l, TheMask)
#RSTwriter(args.outfile_prefix + ".R2l.nc", "R2l", R2l, TheMask)
RSTwriter(args.outfile_prefix + ".R3l.nc", "R3l", R3l, TheMask)
print("done")
netcdf4.write_3d_file(R3l, "R3l", "R3l_yyyy0630-00:00:00.nc", TheMask, compression=True)

