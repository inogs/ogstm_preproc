from bitsea.utilities.argparse_types import existing_dir_path
import argparse
from datetime import datetime

def argument():
    parser = argparse.ArgumentParser(description = 'Executes gzip in parallel')
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = 'The directory containing uncompressed files')

    parser.add_argument(   '--outputdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = 'The directory where you want to dump compressed files')
    parser.add_argument(   '--averagedir', '-a',
                                type = existing_dir_path,
                                required = True,
                                help = 'The directory where you want to dump averaged files')    
    return parser.parse_args()

args = argument()

import glob
import os
import xarray as xr

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


INPUTDIR = args.inputdir
OUTPUTDIR = args.outputdir


filelist=[f for f in INPUTDIR.glob("*-12:00:00.nc") ]
filelist.sort()

coordinates_without_fillvalue = ["nav_lat", "nav_lon",
                "deptht", "deptht_bounds",
                "time_instant", "time_instant_bounds",
                "time_counter", "time_counter_bounds", 
                "time_centered", "time_centered_bounds"]

dateformat_in="%Y-%m-%dT%H:%M:%S.000000000"
dateformat_out="%Y%m%d-%H:%M:%S"

for filename in filelist[rank::nranks]:
    with xr.open_dataset(filename) as ds_file:
        nparts = len(ds_file.time_counter)
        input_var = os.path.basename(filename)[0]  # Estrae la prima lettera (T, U, V, W)
        for it in range(nparts):
            slice_name = str(ds_file.isel(time_counter=it).time_counter.values)
            dt = datetime.strptime(slice_name,dateformat_in)
            outputfile = OUTPUTDIR / f"{input_var}{dt.strftime(dateformat_out)}.nc"
            print("rank %d generates %s" % (rank, outputfile), flush=True)

            ds_slice = ds_file.isel(time_counter=slice(it, it+1))
            for var in ds_slice.variables:
                if var in coordinates_without_fillvalue:
                    ds_slice[var].encoding['_FillValue'] = None
            # for var in ds_slice.variables:
            #     print(f"Variable: {var}, Encoding: {ds_slice[var].encoding['_FillValue'] if '_FillValue' in ds_slice[var].encoding else 'None'  }")
            ds_slice.to_netcdf(outputfile, mode="w", format="NETCDF4")
            ds_slice.close()
        ds_file.mean(dim="time_counter").to_netcdf(OUTPUTDIR / f"{input_var}mean.nc", mode="w", format="NETCDF4")
