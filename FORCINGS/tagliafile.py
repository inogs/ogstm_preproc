import argparse

def argument():
    parser = argparse.ArgumentParser(description = 'Executes gzip in parallel')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'The directory containing uncompressed files')

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = 'The directory where you want to dump compressed files')
    return parser.parse_args()

args = argument()

import glob
import os
from commons.utils import addsep
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False


INPUTDIR =addsep(args.inputdir) 
OUTPUTDIR=addsep(args.outputdir)

filelist=glob.glob(INPUTDIR + "*nc")
filelist.sort()

for filename in filelist[rank::nranks]:
    command="./cut_singlefile_op.sh -i %s -o %s" %( filename, OUTPUTDIR)
    print("rank %d executes %s" % (rank, command), flush=True)
    os.system(command)
    

