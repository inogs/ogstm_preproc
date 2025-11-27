import xarray as xr
import numpy as np
import argparse
import os
import sys

# Parsing command-line arguments
parser = argparse.ArgumentParser(description="Adjust coordinates in NetCDF files.")
parser.add_argument('-i', '--input',    type=str, required=True, help="Input NetCDF file")
parser.add_argument('-o', '--output',   type=str, required=True, help="Output NetCDF file")
parser.add_argument('-v', '--variable', type=str, required=True, help="variable NetCDF file")
args = parser.parse_args()

# File di input e output
file_nc = args.input
OUTDIR = args.output
VAR=args.variable


if VAR.startswith('ch'):
    varmod='P_l'
    OUTDIR='CHL'+OUTDIR
elif VAR.startswith('no'):
    varmod='N3n'
    OUTDIR='NO3'+OUTDIR
elif VAR.startswith('o2'):
    varmod='O2o'
    OUTDIR='O2'+OUTDIR
else:
    sys.exit('new vars must be implemented')


if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

ds = xr.open_dataset(file_nc)
fill_value = np.nan
meshmask_file = "/g100_work/OGS_prodC/OPA/Interim-dev/etc/static-data/MED24_125/meshmask.nc"
meshmask = xr.open_dataset(meshmask_file)
#tmask = meshmask["tmask"].isel(time=0, z=0)
lons = meshmask["nav_lon"][:]
lats = meshmask["nav_lat"][:]
depth   = meshmask["nav_lev"][:] 
new_longitude = lons.values[0,:]


var = ds.variables[VAR][:]

for month in range (0,var.shape[0]):
    MM=month+1
    if len(str(MM)) <=1:
        MM="0"+str( MM) 
    else:
        MM=str(MM)
    tmp = var[month,:,:,:]
    new_shape = list(tmp.shape) 
    new_shape[-1] = 1085
    new_data = np.full(new_shape, fill_value, dtype=tmp.dtype)
    new_data[..., 80:] = tmp.values
    new_ds = xr.Dataset(
            {varmod: (["depth", "latitude", "longitude"], new_data)},
        coords={
            "longitude": (["longitude"], new_longitude),
            "latitude": (["latitude"], lats.values[:,0] ),
            "depth": (["depth"], depth.values),
        }
    )

    new_ds[varmod].attrs[varmod+':_FillValue'] = 1.e+20  # Setting _FillValue
    new_ds[varmod].attrs[varmod+':fillValue'] = 1.e+20   # Setting fillValue
    fileout = os.path.join(OUTDIR, f"ave.yyyy{MM}15-00:00:00.{varmod}.nc")  # Salva nella directory OUTDIR
    new_ds.to_netcdf(fileout,encoding={varmod:{ "zlib":True,  "complevel": 9, "least_significant_digit": 2 }})


"""
# Crea il nuovo dataset con le variabili aggiornate
ds_new = ds.drop_vars("longitude")
ds_new = ds_new.assign(updated_vars)

# Salva il file modificato
ds_new.to_netcdf(output_nc)


meshmask_file = "/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022_v1/wrkdir/MODEL/meshmask.nc"
meshmask = xr.open_dataset(meshmask_file)
tmask = meshmask["tmask"].isel(time=0, z=0)
lons = meshmask["nav_lon"][:]
lats = meshmask["nav_lat"][:]
meshmask.close()
"""
