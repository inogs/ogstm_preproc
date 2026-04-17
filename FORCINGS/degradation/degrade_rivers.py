import numpy as np
#from sys import exit
import xarray as xr
from glob import glob
from pathlib import Path
from argparse import ArgumentParser

from itertools import chain
# my stuff
import degrade_mesh as dm
from commons import load_parameters
from degrade_forcings import load_mesh_light, vwmean2d
from bitsea.commons.mask import Mask
import forcingswriter as FW

'''
degrades resolution of river load files TIN_*.nc
this is because /ogstm_preproc/BC/EFAS/rivers.xml
has ij of destination grid hardcoded, so it only works on the 1/24 deg grid.
Here must take care the degradation algorithm for loads here
is consistent with that used to degrade runoff ("sorunoff") 
in degrade forcings.py
'''

def argument():
    parser = ArgumentParser()
    parser.add_argument('--yamlfile','-y',
                        type=str,
                        required=True,
                        help='file with all parameters')
    return parser.parse_args()

def get_flist_riv(Params:dict):
    rivdir = Params['rivdir']
    y0 = Params['y0']
    yE = Params['yE']
    flist = [glob(f'{rivdir}/TIN_{YYYY}*.nc') for YYYY in range(y0, yE+1)]
    flist = sorted(list(chain.from_iterable(flist)))
    return [Path(f) for f in flist]

def get_area(maskfile):
    '''
    givern a meshmask, returns cell surface area (2D) on W (T) grid
    '''
    M = xr.open_dataset(maskfile)
    e1t = M['e1t'].values[:]
    e2t = M['e2t'].values[:]
    At = xr.DataArray(e1t * e2t, name='At', dims=('time','z_a','y','x'))
    #tmask = M['tmask'].values[:]
    #tmask_0 = tmask[0,0,:,:].astype(bool)
    M.close()
    return At

def load_rfile(rfile, ndeg=1):
    '''
    all 2D vars, no coordinates!
    dims: (lat, lon)
    '''
    F = xr.open_dataset(rfile)
    F = F.rename({'lat':'y', 'lon':'x'})
    F1 = {}
    #
    for vv in list(F.variables.keys()):
        X = F[vv]
        X = X.where(X != -1.0, 0.0) #take away the -1, else they get in the average
        F1[vv] = dm.xpnd_wrap(X, 'edge', ndeg)
    #
    F1 = xr.Dataset(F1)
    #F1 = F1.rename({'y':'lat', 'x':'lon'}) #probably not needed!
    return F1

def degrade_RIV(R, tmask_in, Warea, Mask_out, outfile, ndeg=1):
    '''
    degrades resolution of T-grid RIV file (2D)
    R        : xarray Dataset, t-grid with river loads
    tmask_in : xarray DataArray, reshaped in blocks tmask of original mesh
    Warea    : xarray DataArray, reshaped in blocks cell surface area on t-grid
    Mask_out : Mask object, degraded mesh mask
    outfile  : str, path to output file
    ndeg     : int, degradation factor
    '''
    Rd = {} 
    for vv in list(R.variables.keys()):
        X = dm.reshape_blocks(R[vv], ndeg)
        Xd = vwmean2d(X, tmask_in, Warea, Mask_out, mask_weight=True).squeeze()
        Xd = Xd.where(Xd != 0.0, -1.0) #put the -1s back at sea
        Rd[vv] = Xd
    Rd = xr.Dataset(Rd)
    return Rd

if __name__=='__main__':
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nranks = comm.size
    except:
        comm = None
        rank = 0
        nranks = 1
    #
    #yamlfile = 'degrade_rivers.yaml'; Params = load_parameters(yamlfile)
    Params = load_parameters(argument().yamlfile)
    #
    maskfile = Path(Params['maskfile'])
    maskfile_d = Path(Params['maskfile_d'])
    '''
    WARNING! THESE ARE THE MESHES WITH ATLANTIC BUFFER,
    WHERE "sorunoff" LIVES, while the TIN_*.nc files
    LIVE ON THE OTHER MESH, SO THEY MUST MATCH!!!
    goddammit!!!
    '''
    rivdir =  Path(Params['rivdir'])
    outdir = Path(Params['outdir'])
    outdir.mkdir(exist_ok=True, parents=True)
    y0 = Params['y0']
    yE = Params['yE']
    ndeg = Params['ndeg']
    #
    M = load_mesh_light(maskfile, ndeg)
    tmask_in = dm.reshape_blocks(M['tmask'].astype(bool), ndeg)
    Mask_out = Mask.from_file(maskfile_d)
    At = get_area(maskfile)
    At = dm.xpnd_wrap(At, 'interp', ndeg)
    Warea = dm.reshape_blocks(At, ndeg)
    #
    flist = get_flist_riv(Params)
    #
    for rfile in flist[rank::nranks]:
        '''
        HERE MUST PUT SOMETHING SO THAT RUNOFF AND LOAD LINE UP!
        '''
        outfile = outdir / rfile.name
        R = load_rfile(rfile, ndeg)
        Rd = degrade_RIV(R, tmask_in, Warea, Mask_out, outfile, ndeg)
        Rd.to_netcdf(outfile, format='NETCDF4_CLASSIC')
