import argparse
from bitsea.commons.mask import Mask
from IC import RSTwriter
from bitsea.commons.dataextractor import DataExtractor

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates restarts for R8 vars, reading from R6 ones
    Applied rule is : R8 = R6 * 0.1
    ''')
    parser.add_argument(   '--file_prefix', '-o',
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

TheMask = Mask.from_file(args.maskfile)

for var in ["c","n","p","s"]:
    inputfile = f"{args.file_prefix}.R6{var}.nc"
    outfile = f"{args.file_prefix}.R8{var}.nc"
    R6 = DataExtractor(TheMask,inputfile,f"TRNR6{var}").values
    R8 = R6 * 0.1
    RSTwriter(outfile, f"R8{var}", R8,TheMask)

