# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES

import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates a restart for O5c, with three values: top, middle and bottom
    ''')
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = False,
                                default = "RST.19950101-00:00:00.O5c.nc",
                                help = 'Output directory of generated restarts'
                                )
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )
    return parser.parse_args()

args = argument()



import numpy as np
from commons.mask import Mask
from IC import RSTwriter


TheMask=Mask(args.maskfile)


jpk, jpj, jpi=TheMask.shape

time=1

nav_lev=TheMask.zlevels

#    double PIC(time, z, y, x) ;
PIC = np.ones((jpk,jpj,jpi),np.double)*1.e+20;
for jk in range(jpk):
    for jj in range(jpj):
        for ji in range(jpi):
            if TheMask.mask[jk,jj,ji]:
                if nav_lev[jk] < 150 :
                    PIC[jk,jj,ji]  = 2.0
                elif nav_lev[jk] >= 150 and nav_lev[jk] < 300:
                    PIC[jk,jj,ji]  = 1.0
                else:
                    PIC[jk,jj,ji]  = 1.e-11

##############################################################
# write meshmask netcdf file !
##############################################################
var="O5c"

RSTwriter(args.outfile, var, PIC, TheMask)


