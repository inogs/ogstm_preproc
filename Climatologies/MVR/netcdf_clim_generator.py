import xarray as xr
import numpy as np
import argparse
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path


# Parsing command-line arguments
parser = argparse.ArgumentParser(description="Adjust coordinates in NetCDF files.")
parser.add_argument('-i', '--input',    type=existing_file_path, required=True, help="Input NetCDF file")
parser.add_argument('-o', '--outdir',   type=existing_dir_path, required=True, help="string for output directory")
parser.add_argument('-v', '--variable', type=str, required=True, choices=['chl_avg', 'no3_avg', 'o2_avg'], help="variable NetCDF file")
parser.add_argument('-m', '--maskfile', type=existing_file_path, required=True, help="meshmask file")
parser.add_argument('-d', '--least_significant_digit', type=int, required=False, default=0, help="least significant digit for compression")
args = parser.parse_args()

VAR=args.variable

VARMOD={
    'chl_avg':'P_l',
    'no3_avg':'N3n',
    'o2_avg':'O2o'
}
varmod=VARMOD[VAR]
if args.least_significant_digit == 0:
    least_significant_digit = None
else:
    least_significant_digit = args.least_significant_digit

TheMask = Mask.from_file(args.maskfile)
jpk,jpi,jpj = TheMask.shape

with xr.open_dataset(args.input) as ds:
    var = ds.variables[VAR][:]
fill_value = 1.e+20


for month in range(12):
    fileout = args.outdir / f"ave.yyyy{month+1:02d}15-00:00:00.{varmod}.nc"
    print(f"Saving to {fileout} with least_significant_digit = {least_significant_digit}", flush=True)
    cmems_clim = var[month,:,:,:]

    new_data = np.full(TheMask.shape, fill_value, dtype=np.float32)
    new_data[:,:,80:] = cmems_clim.values

    new_ds = xr.Dataset(
            {varmod: (["depth", "latitude", "longitude"], new_data)},
        coords={
            "longitude": (["longitude"], TheMask.lon),
            "latitude": (["latitude"], TheMask.lat),
            "depth": (["depth"], TheMask.zlevels),
        }
    )

    new_ds.to_netcdf(fileout,
                     encoding={varmod:{ 
                         "zlib":True,  
                         "_FillValue": fill_value, 
                         "complevel": 9, 
                         "least_significant_digit": least_significant_digit,
                         }
                         }
                     )
    new_ds.close()