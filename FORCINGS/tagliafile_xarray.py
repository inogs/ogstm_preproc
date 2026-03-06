from bitsea.utilities.argparse_types import existing_dir_path
import argparse
from datetime import datetime

def argument():
    parser = argparse.ArgumentParser(description = 'Executes gzip in parallel')
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = 'Directory with original CMCC files with 4 time frames per day')

    parser.add_argument(   '--outputdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = 'Directory with cut files with 1 time frame per file, every 6 hours')
    parser.add_argument(   '--averagedir', '-a',
                                type = existing_dir_path,
                                required = False,
                            help = """Directory with 24h averaged files. 
                                      If not specified, averaged files will not be generated
                                      """)    
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
AVERAGEDIR = args.averagedir

filelist=[f for f in INPUTDIR.glob("*-12:00:00.nc") ]
filelist.sort()

coordinates_without_fillvalue = ("nav_lat", "nav_lon",
                "deptht", "deptht_bounds",
                "time_instant", "time_instant_bounds",
                "time_counter", "time_counter_bounds", 
                "time_centered", "time_centered_bounds")
if AVERAGEDIR is not None:
    valid_nc4_keys = ('dtype', 'zlib', 'complevel', 'fletcher32', 'contiguous',
                        'chunksizes', '_FillValue', 'shuffle', 'endian',
                        'least_significant_digit', 'significant_digits',
                        'compression', 'quantize_mode', 'szip_coding',
                        'szip_pixels_per_block', 'blosc_shuffle', 'missing_value')

dateformat_in="%Y-%m-%dT%H:%M:%S.000000000"
dateformat_out="%Y%m%d-%H:%M:%S"

for filename in filelist[rank::nranks]:
    with xr.open_dataset(filename) as ds_file:
        nparts = len(ds_file.time_counter)
        input_var = os.path.basename(filename)[0]  # Estrae la prima lettera (T, U, V, W)
        for it in range(nparts):
            # genera il nome del file di output che indica l'ora centrale della finestra temporale
            slice_name = str(ds_file.isel(time_counter=it).time_counter.values)
            dt = datetime.strptime(slice_name,dateformat_in)
            outputfile = OUTPUTDIR / f"{input_var}{dt.strftime(dateformat_out)}.nc"
            print("rank %d generates %s" % (rank, outputfile), flush=True)
            # estrae la fetta temporale corrispondente
            ds_slice = ds_file.isel(time_counter=slice(it, it+1))
            for var in ds_slice.variables:
                if var in coordinates_without_fillvalue:
                    ds_slice[var].encoding['_FillValue'] = None
            # salva il file
            ds_slice.to_netcdf(outputfile, mode="w", format="NETCDF4")
            ds_slice.close()
        if AVERAGEDIR is not None:
            outputfile = AVERAGEDIR / os.path.basename(filename)
            print("rank %d generates %s" % (rank, outputfile), flush=True)
            # Calcola la media temporale su un arco temporale di 24h e salva il file
            ds_mean = ds_file.mean(dim="time_counter", keep_attrs=True).expand_dims("time_counter")

            # ---- settings to get correctly masked values in ncdump -h
            for var in ds_mean.variables:
                if var.startswith(('vo', 'so')):
                    ds_mean[var].encoding['_FillValue'] = 1.e20
                    ds_mean[var].encoding['missing_value'] = 1.e20
                if var.startswith('depth'):
                    ds_mean[var].encoding['_FillValue'] = None
            # Preserve each variable's original encoding, then enable compression and raise complevel to 4 for data variables
            encoding = {}
            for var in ds_mean.variables:
                enc = {k: v for k, v in ds_mean[var].encoding.items() if k in valid_nc4_keys}
                if var in ds_mean.data_vars:
                    enc['zlib'] = True
                    enc['complevel'] = 5
                    enc['shuffle'] = True
                encoding[var] = enc
            # -------------------------------------------------
            ds_mean.to_netcdf(outputfile, mode="w", format="NETCDF4", encoding=encoding)
            ds_mean.close()