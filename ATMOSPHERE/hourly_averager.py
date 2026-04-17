from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_averagers import TimeAverager2D
from bitsea.commons.mask import Mask
from bitsea.commons import netcdf4
from pathlib import Path
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path

import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Interpolates with nearest algorithm restarts between very similar meshes
    Meshes are supposed to have the same definition and be different from some land/sea points
 
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = 'original restarts')

    parser.add_argument(   '--outdir', '-o',
                                type = Path,
                                required = True,
                                help = 'Output directory of generated restarts'
                                )

    parser.add_argument(   '--maskfile','-m',
                                type = existing_file_path,
                                required = False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V4.1/meshmask.nc",
                                help = 'Path of the mask file'
                                )    

    return parser.parse_args()

args = argument()

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

TheMask = Mask.from_file(args.maskfile)
DIR = Path(args.inputdir)
OUTDIR= Path(args.outdir)

OUTDIR.mkdir(parents=True, exist_ok=True)


TL = TimeList.fromfilenames(None, DIR, "??/atm*", prefix="atm.")
Hourly_reqs= TL.getHourlyList(hours=6,start_hour=3)

varlist=["sp","msl","t2m","wsp10","d2m","tcc","tco3","tclw"]

for hourly_req in Hourly_reqs[rank::nranks]:
    ii,w = TL.select(hourly_req)
    filelist=[TL.filelist[k] for k in ii]
    yyyy=hourly_req.centertime.strftime('%Y')
    mm =hourly_req.centertime.strftime('%m')
    outdir = OUTDIR / yyyy / mm
    (OUTDIR / yyyy).mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile= outdir / hourly_req.centertime.strftime("atm.%Y%m%d-%H:%M:%S.nc")

    print(f"Rank {rank} generates {outfile}",flush=True)
    for var in varlist:
        if var == "wsp10":
            U=TimeAverager2D(filelist,w,"u10", TheMask)
            V=TimeAverager2D(filelist,w,"v10", TheMask)
            M2d=(U**2+V**2)**0.5
        else:
            M2d=TimeAverager2D(filelist,w,var, TheMask)
        netcdf4.write_2d_file(M2d,var,outfile,TheMask)


