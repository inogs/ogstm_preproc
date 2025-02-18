# this script requires netcdf 4
# to load them it is needed the following virtual environment:
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh
# LOAD PACKAGES

import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates a restart for R1l, defined as R3c*0.02
    ''')
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = False,
                                default = "RST.19950101-00:00:00.R1l.nc",
                                help = 'Output directory of generated restarts'
                                )
    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = 'Path of a R3c RST file'
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
from bitsea.commons.dataextractor import DataExtractor

TheMask=Mask.from_file(args.maskfile)


R3c = DataExtractor(TheMask, args.inputfile, "TRNR3c").values
R1l = R3c*0.02


RSTwriter(args.outfile, "R1l", R1l, TheMask)


