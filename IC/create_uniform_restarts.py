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
from IC import RSTwriter
import numpy as np


TheMask=Mask(args.maskfile)


jpk,jpj,jpi=TheMask.shape
R1l = np.ones((jpk,jpj,jpi), np.float64)*0.0
R2l = np.ones((jpk,jpj,jpi), np.float64)*0.0
R3l = np.ones((jpk,jpj,jpi), np.float64)*1.0
R3c = np.ones((jpk,jpj,jpi), np.float64)*600.0

RSTwriter(args.outfile_prefix + ".R1l.nc", "R1l", R1l, TheMask)
RSTwriter(args.outfile_prefix + ".R2l.nc", "R2l", R2l, TheMask)
RSTwriter(args.outfile_prefix + ".R3l.nc", "R3l", R3l, TheMask)
RSTwriter(args.outfile_prefix + ".R3c.nc", "R3c", R3c, TheMask)

