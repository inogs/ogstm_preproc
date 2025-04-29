import numpy as np
from bitsea.commons import netcdf4
from bitsea.commons.mask import Mask

def generate(mask,filename):
    val = 1./3600
    jpk,jpj,jpi = mask.shape
    V = np.ones((jpk,jpj,jpi),np.float32)*val
    V[~mask.mask] = 1.e+20
    netcdf4.write_3d_file(V, "reR3l", filename, mask, compression=True)

if __name__=="__main__":
    import argparse
    def argument():
        parser = argparse.ArgumentParser(description = '''
        Creates R3l restoration file                    
                                   ''')
        parser.add_argument(   '--outfile', '-o',
                                    type = str,
                                    required = True,
                                    help="bounmask.nc"
                                    )
        parser.add_argument(   '--maskfile', '-m',
                                    type = str,
                                    required = True,
                                    help ='/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc'
                                    )
        return parser.parse_args()
    args = argument()
   
    TheMask = Mask.from_file(args.maskfile)
    generate(TheMask,args.outfile)
