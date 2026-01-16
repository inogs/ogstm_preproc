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
    Generates restart and nudging file for R3l
    ''')
    parser.add_argument(   '--rst_file', '-r',
                                type = str,
                                required = True,
                                default = "RST.19950101-00:00:00.R3l.nc",
                                help = "Path of restart file"
                                )
    parser.add_argument(   '--ndg_file', '-n',
                                type = str,
                                required = True,
                                default = "R3l.yyyy0630-00:00:00.R3l.nc",
                                help = 'Path of nudging file'
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


RSTwriter(args.rst_file, "R3l", R3l_s, TheMask)
netcdf4.write_3d_file(R3l_s, "R3l", args.ndg_file, TheMask, compression=True)

