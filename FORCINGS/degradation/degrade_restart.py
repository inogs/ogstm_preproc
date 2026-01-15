import numpy as np
import xarray as xr
import netCDF4
from argparse import ArgumentParser
from glob import glob
#
import degrade_mesh as dm
from commons import degrade_wrap, load_parameters
import regridding as rg
from pathlib import Path
from bitsea.commons.mask import Mask
#from IC import RSTwriter
#from previous.restart_interpolator import RSTwriter

'''
degrades the horizontal resolution of BFM restarts by an integer value
joins boxes of ndeg x ndeg cells, does volume weighted mean.
needs both original and degraded-resolution meshmask 
(run degrade mesh.py before)

usage: 
python degrade_mesh.py -y degrade_mesh.yaml

yaml file example:
mesh_in: '/path/to/meshmask_original.nc'
mesh_out: '/path/to/meshmask_degraded.nc'
indir: '/path/to/restarts/rst.nc' #also with wildcards '*'
outdir: '/path/to/output_dir/'
infile: 'RST.*-00:00:00.__VNAME__.nc' #__VNAME__ gets replaced and is needed to read the variable
ndeg: 6 #size of boxes for resolution degradation
'''

def argument():
    parser = ArgumentParser()
    parser.add_argument('--yamlfile','-y',
                        type=str,
                        required=True,
                        help='file with all parameters')
    return parser.parse_args()


def get_Volume(M):
    '''
    gets cell volume for weighted mean
    '''
    V0 = M['e1t'].values[:] * M['e2t'].values[:] * M['e3t_0']
    nanmask = M['tmask'].values[:]
    nanmask[nanmask==0.0] = np.nan
    V0 = V0 * nanmask #else you see border effects at the coast
    return V0


def load_xpnd_rst(infile, vname, mask:Mask, ndeg=1):
    '''
    loads  the data array from a variable in the restart and expands it in 'edge' mode
    Arguments:
        infile : path to restart file
        vname : bfm variable name
        mask : Mask object of the fine mesh

    Returns:
    RST : xarray DataArray, with 3D variable already expanded
         in 'edge' mode, nan 
    '''
    vname = "TRN"+vname
    with xr.open_dataset(infile) as D:
        A = D[vname]
        A.values[0,:][~mask.mask] = np.nan
        RST = dm.xpnd_wrap(A, 'edge', ndeg)
    return RST


def degrade_bgc(DI, V0, Maskout:Mask, ndeg=1):
    '''
    degrades resolution of BFM restart file
    Arguments:
    D1 : xarray DataArray, with 3D variable already expanded
         in 'edge' mode
    V0: array DataArray of the Volume of fine mesh
    Maskout : mask of coarse mesh

    Returns:
    R : xarray Dataset, degraded variable

    '''
    
    A = degrade_wrap(DI, dm.vwmean, V0, ndeg)

    return A

def get_flist(Params):
    indir = Params['indir']
    infile = Params['infile']
    flist = glob(indir + infile.replace('__VNAME__', '*'))
    itrbl = [(ff.split('.')[-2], Path(ff)) for ff in flist]
    return itrbl

def RSTwriter(outfile, var, rst, TheMask):
    rst[~TheMask.mask] = 1.e+20
    jpk, jpj, jpi = TheMask.shape
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',1)

    TRN   = 'TRN' + var;
    ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = TheMask.ylevels
    ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = TheMask.zlevels
    ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
    ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x') ); ncvar[:] = rst;

    setattr(ncOUT.variables[TRN]   ,'missing_value',1.e+20                             );
    setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
    setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
    ncOUT.close()
    return


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
    Params = load_parameters(argument().yamlfile)
    mesh_in = Path(Params['mesh_in'])
    mesh_out = Path(Params['mesh_out'])
    outdir = Path(Params['outdir'])
    ndeg = Params['ndeg']

    M = dm.load_mesh(mesh_in, ndeg)
    V0 = get_Volume(M)

    itrbl = get_flist(Params)

    Maskout = Mask.from_file(mesh_out)
    Mask_in = Mask.from_file(mesh_in)

    for vname, fname in itrbl[rank::nranks]:
        outfile = outdir / fname.name
        print(f'degrado: {vname}')
        # 
        DI = load_xpnd_rst(fname, vname, Mask_in,ndeg)
        Dd = degrade_bgc(DI, V0, Maskout, ndeg)

        print("eccolo")
        v = Dd.values[0,:]
        #v[~Maskout.mask] = 1.0e20
        waterpoints=v[Maskout.mask]
        print(f' post-degradation min/max/nans/ {waterpoints.min()}, {waterpoints.max()}', np.sum(np.isnan(waterpoints)))
        RSTwriter(outfile, vname, Dd.values[0,:], Maskout)

