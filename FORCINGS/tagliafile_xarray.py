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

coordinates_without_fillvalue = ["nav_lat", "nav_lon",
                "deptht", "deptht_bounds",
                "time_instant", "time_instant_bounds",
                "time_counter", "time_counter_bounds", 
                "time_centered", "time_centered_bounds"]
for filename in filelist[rank::nranks]:
    with xr.open_dataset(filename) as ds_file:
        print("rank %d executes cut of %s" % (rank, filename), flush=True)
        for it in range(4):
            slice_name = str(ds_file.isel(time_counter=it).time_counter.values)
            slice_name = slice_name.split('.')[0]
            slice_time = slice_name.split('T')[1]
            slice_date = slice_name.split('T')[0].replace('-', '')
            outputfile = OUTPUTDIR + 'T' + slice_date + "-" + slice_time + ".nc"
            ds_slice = ds_file.isel(time_counter=slice(it, it+1))
            for var in ds_slice.variables:
                if var in coordinates_without_fillvalue:
                    ds_slice[var].encoding['_FillValue'] = None
            # for var in ds_slice.variables:
            #     print(f"Variable: {var}, Encoding: {ds_slice[var].encoding['_FillValue'] if '_FillValue' in ds_slice[var].encoding else 'None'  }")
            ds_slice.to_netcdf(outputfile, mode="w", format="NETCDF4")
            ds_slice.close()
            print("rank %d generates %s" % (rank, outputfile), flush=True)


    
    

