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
import xarray as xr
from bitsea.commons.utils import addsep
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
    ds_file = xr.open_dataset(filename)
    print("rank %d executes cut of %s" % (rank, filename), flush=True)
    for it in range(4):
        slice_name = str(ds_file.isel(time_counter=it).time_counter.values)
        slice_name = slice_name.split('.')[0]
        slice_time = slice_name.split('T')[1]
        slice_date = slice_name.split('T')[0].replace('-', '')
        outputfile = OUTPUTDIR + 'T' + slice_date + "-" + slice_time + ".nc"
        ds_file.isel(time_counter=it).to_netcdf(outputfile, mode="w", format="NETCDF4")
        # ds_file.isel(time_counter=it).to_netcdf(outputfile, mode="w", format="NETCDF4")
        print("rank %d generates %s" % (rank, outputfile), flush=True)
    ds_file.close()

    
    

