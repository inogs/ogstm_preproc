import argparse
from datetime import datetime

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

# Crea la directory output se non esiste (solo rank 0)
if rank == 0 and not os.path.exists(OUTPUTDIR.rstrip('/')):
    os.makedirs(OUTPUTDIR.rstrip('/'), exist_ok=True)
    print(f"Created output directory: {OUTPUTDIR}", flush=True)

if isParallel:
    comm.Barrier()  # Sincronizza prima di iniziare

filelist=glob.glob(INPUTDIR + "*-12:00:00.nc")
filelist.sort()

coordinates_without_fillvalue = ["nav_lat", "nav_lon",
                "deptht", "deptht_bounds",
                "time_instant", "time_instant_bounds",
                "time_counter", "time_counter_bounds", 
                "time_centered", "time_centered_bounds"]

dt_format_in="%Y-%m-%dT%H:%M:%S.000000000"
dt_format_out="%Y%m%d-%H:%M:%S"

for filename in filelist[rank::nranks]:
    with xr.open_dataset(filename) as ds_file:
        nparts = len(ds_file.time_counter)
        input_var = os.path.basename(filename)[0]  # Estrae la prima lettera (T, U, V, W)
        print("rank %d executes cut of %s in %d parts" % (rank, filename, nparts), flush=True)
        for it in range(nparts):
            slice_name = str(ds_file.isel(time_counter=it).time_counter.values)
            dt = datetime.strptime(slice_name,dt_format_in) # format check
            print(f"rank {rank} processing slice: {slice_name}", flush=True)

            dt_name = dt.strftime(dt_format_out)
            outputfile = OUTPUTDIR + input_var + dt_name + ".nc"
            print("rank %d generates %s" % (rank, outputfile), flush=True)

            ds_slice = ds_file.isel(time_counter=slice(it, it+1))
            for var in ds_slice.variables:
                if var in coordinates_without_fillvalue:
                    ds_slice[var].encoding['_FillValue'] = None
            # for var in ds_slice.variables:
            #     print(f"Variable: {var}, Encoding: {ds_slice[var].encoding['_FillValue'] if '_FillValue' in ds_slice[var].encoding else 'None'  }")
            ds_slice.to_netcdf(outputfile, mode="w", format="NETCDF4")
            ds_slice.close()
            print("rank %d generates %s" % (rank, outputfile), flush=True)
