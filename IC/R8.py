import argparse
from bitsea.utilities.argparse_types import existing_file_path, existing_dir_path
from bitsea.commons.mask import Mask
from IC import RSTwriter
from bitsea.commons.dataextractor import DataExtractor

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates restarts for R8 vars, reading from R6 ones
    Applied rule is : R8 = R6 * 0.1
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = 'Path of the input directory'
                                )
    parser.add_argument(   '--date', '-d',
                                type = str,
                                required = True,
                                help = 'Date of the restarts'
                                )    
    parser.add_argument(   '--maskfile','-m',
                                type = existing_file_path,
                                required = True,
                                help = 'Path of the mask file'
                                )    
    return parser.parse_args()

args = argument()

TheMask = Mask.from_file(args.maskfile)

for var in ["c","n","p","s"]:
    inputfile = args.inputdir / f"RST.{args.date}.R6{var}.nc"
    outfile = args.inputdir / f"RST.{args.date}.R8{var}.nc"
    print(f"Processing {inputfile} to generate {outfile}")
    R6 = DataExtractor(TheMask,inputfile,f"TRNR6{var}").values
    R8 = R6 * 0.1
    RSTwriter(outfile, f"R8{var}", R8,TheMask)

