# to load them it is needed the following virtual environment:
# source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh
# usage:
# opfx='/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/SCRIPTS/RST.19950101-00:00:00'
# maskfile='/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022_v6/wrkdir/MODEL/meshmask.nc'
# csvfile='/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/SCRIPTS/mean_R3l_from_RRS412.csv'
# python create_N-basins_restarts.py -o $opfx -m $maskfile -c $csvfile

from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.basins import V2
from bitsea.basins.basin import ComposedBasin
from bitsea.commons.interpolators import shift
from IC import RSTwriter
import numpy as np
import pandas as pd
from bitsea.commons import netcdf4
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates a restart for R1l, R2l R3l, R3c
    ''')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = "out/dir/",
                                help = 'Path of outfiles'
                                )
    parser.add_argument(   '--rst_file', '-r',
                                type = str,
                                required = True,
                                default = "/RST.19950101-00:00:00",
                                help = 'outfiles without .R1l.nc'
                                )
    parser.add_argument(   '--ndg_file', '-n',
                                type = str,
                                required = True,
                                default = "",
                                help = 'outfiles without .R1l.nc'
                                )
    parser.add_argument(   '--maskfile','-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )   
    parser.add_argument(   '--csvfile','-c',
                                type = str,
                                required = True,
                                help = 'Path of the csv file with surface R3l values'
                                )
    return parser.parse_args()

def CDOM_profile(surface_value, bottom_value, depth, depth0):
    k = 0.05
    res = surface_value + (bottom_value - surface_value) / (1.0 + np.exp(-k * (depth - depth0)))
    return res

def smoother3D(mask, RST):
    #NsmoothingCycles = 20
    NsmoothingCycles = 2000
    jpk,jpj,jpi = mask.shape
    RST[mask ==0 ] = np.nan
    auxs = np.zeros((5,jpk,jpj,jpi),np.float32)
    for ss in range(NsmoothingCycles):
        print(f'{ss}/{NsmoothingCycles}')
        for j in range(jpk):
           auxs[0,j,:,:] = RST[j,:,:]
           auxs[1,j,:,:] = shift(RST[j,:,:],1,'r')
           auxs[2,j,:,:] = shift(RST[j,:,:],1,'l')
           auxs[3,j,:,:] = shift(RST[j,:,:],1,'u')
           auxs[4,j,:,:] = shift(RST[j,:,:],1,'d')
        RST = np.nanmean(auxs,axis=0)
    return RST

def smoother2D(mask,RST):
    NsmoothingCycles = 2000
    jpj,jpi = mask.shape
    RST[mask ==0 ] = np.nan
    auxs = np.zeros((5,jpj,jpi),np.float32)
    for ss in range(NsmoothingCycles):
        print(f'{ss}/{NsmoothingCycles}')
        auxs[0,:,:] = RST
        auxs[1,:,:] = shift(RST,1,'r')
        auxs[2,:,:] = shift(RST,1,'l')
        auxs[3,:,:] = shift(RST,1,'u')
        auxs[4,:,:] = shift(RST,1,'d')
        RST = np.nanmean(auxs,axis=0)
    return RST

args = argument()
TheMask = Mask.from_file(args.maskfile)
tmask_3D = TheMask._data_array.astype(float)
csvfile = args.csvfile

Surf_Values = pd.read_csv(csvfile, index_col=0)

jpk, jpj, jpi = TheMask.shape
R3l_0 = np.ones((1, jpj, jpi), np.float64) * 1.0
R3l = np.ones((jpk, jpj, jpi), np.float64) * 1.0

nav_lev=TheMask.zlevels

for sub in V2.P:
    bname = sub.name
    if bname=='med':
        continue
    print(bname)
    surf_val = Surf_Values.loc[bname].R3l
    curr_mask = SubMask(sub, mask=TheMask)
    R3l_0[0,curr_mask[0,:,:]] = surf_val

print('smoothing...')
R3l_s = smoother2D(tmask_3D[0,:,:], R3l_0[0,:,:])
#R3l_s = R3l_0[0,:,:]

R3l_s = CDOM_profile(R3l_s[np.newaxis,:,:], 1.0, nav_lev[:,np.newaxis,np.newaxis], 150.0)

outdir = args.outdir
outfile_rst = outdir+args.rst_file+'.R3l.nc'
outfile_ndg = outdir+args.ndg_file
RSTwriter(outfile_rst, "R3l", R3l_s, TheMask)
print("done")
netcdf4.write_3d_file(R3l_s, "R3l", outfile_ndg, TheMask, compression=True)

