import numpy as np
import xarray as xr
from argparse import ArgumentParser
import yaml
from glob import glob
#
import degrade_mesh as dm
import degrade_forcings as df

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

def load_parameters():
    yamlfile=argument().yamlfile
    # yamlfile = 'degrade_restart.yaml'
    with open(yamlfile) as f:
        Params = yaml.load(f, Loader=yaml.Loader)
    return(Params)

def get_Volume(M):
    '''
    gets cell volume for weighted mean
    '''
    V0 = M['e1t'].values[:] * M['e2t'].values[:] * M['e3t_0']
    nanmask = M['tmask'].values[:]
    nanmask[nanmask==0.0] = np.nan
    V0 = V0 * nanmask #else you see border effects at the coast
    return V0

def load_rst(infile, vname, ndeg=1):
    '''
    loads the data array from a variable in the restart
    '''
    D = xr.open_dataset(infile)
    D1 = {}
    #
    D1[vname] = dm.xpnd_wrap(D[vname], 'edge', ndeg)
    #
    D1['nav_lat'] = dm.xpnd_wrap(D['nav_lat'], 'interp', ndeg)
    D1['nav_lon'] = dm.xpnd_wrap(D['nav_lon'], 'interp', ndeg)
    #D1[vname].coords['nav_lat'].values[:] = D1['nav_lat'].values[:]
    #D1[vname].coords['nav_lon'].values[:] = D1['nav_lon'].values[:]
    D1['time'] = D['time']
    D1['nav_lev'] = D['nav_lev']
    #
    D1 = xr.Dataset(D1)
    D1 = D1.assign_attrs(D.attrs)
    D.close()
    return D1

def init_rst(C):
    '''
    lat, lon are the same for all files
    '''
    R = {}
    R['nav_lon'] = C['glamt']
    R['nav_lat'] = C['gphit']
    return R

def degrade_bgc(DI, V0, C, vname, ndeg=1):
    '''
    degrades resolution of BFM restart file
    '''
    R = init_rst(C)
    R[vname] = df.degrade_wrap(DI[vname], dm.vwmean, V0, ndeg)
    R['nav_lev'] = DI['nav_lev']
    R = xr.Dataset(R)
    R = R.assign_attrs(DI.attrs)
    return R

def get_flist(Params):
    indir = Params['indir']
    infile = Params['infile']
    flist = glob(indir+infile.replace('__VNAME__', '*'))
    itrbl = [(ff.split('.')[-2], ff) for ff in flist]
    return itrbl

if __name__=='__main__':
    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nranks = comm.size
    except:
        comm = None
        rank = 0
        nranks = 1
    #
    Params = load_parameters()
    mesh_in = Params['mesh_in']
    mesh_out = Params['mesh_out']
    indir = Params['indir']
    infile = Params['infile']
    outdir = Params['outdir']
    ndeg = Params['ndeg']
    #
    if rank==0:
        M = df.load_mesh_light(mesh_in, ndeg)
        C = df.load_coords_degraded(mesh_out)
        V0 = get_Volume(M)
    else:
        M = None
        C = None
        V0 = None
    if nranks > 1:
        M = comm.bcast(M, root=0)
        C = comm.bcast(C, root=0)
        V0 = comm.bcast(V0, root=0)
    # HERE LOOP ON FILES AND VARIABLES
    itrbl = get_flist(Params)
    for vname, fname in itrbl[rank::nranks]:
        vname = 'TRN'+vname
        print(f'degrado: {vname}')
        # 
        DI = load_rst(fname, vname, ndeg)
        Dd = degrade_bgc(DI, V0, C, vname, ndeg)
        outfile = outdir+fname.split('/')[-1]
        dm.dump_netcdf(Dd, outfile)
