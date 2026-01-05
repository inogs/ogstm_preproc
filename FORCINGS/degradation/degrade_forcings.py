import numpy as np
import xarray as xr
from glob import glob
import os
from mpi4py import MPI
from argparse import ArgumentParser
import yaml
from itertools import chain
# my stuff
import degrade_mesh as dm
from commons import dump_netcdf 

'''
degrades resolution of OGSTM physics forcings
on T, U, V, W grids
needs both original and degraded resolution meshmask.nc
(run degrade_mesh.py before)
input parameters (like directories etc.) go into a .yaml file
usage: python degrade_forcings.py -y degrade_forcings.yaml

yaml file example:
maskfile: '/path/to/meshmask_in.nc' #original mesh
maskfile_d: '/path/to/coarse/meshmask_out.nc' #degraded mesh
ffdir: '/path/to/wrkdir/MODEL/FORCINGS/'
outdir: '/path/to/output/'
ndeg: 6   #nr of cells to join in i, j
y0: 2000  #start year
yE: 2020  #end year

(takes about 3min per year of daily files on 1 node (256GB) and 16 cores,
if it doesn't die by OOM, which sometimes does and sometimes not. 
If it does go with 8 cores instead)
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
    # yamlfile = 'degrade_forcings.yaml'
    with open(yamlfile) as f:
        Params = yaml.load(f, Loader=yaml.Loader)
    return(Params)

def load_mesh_light(maskfile, ndeg=1):
    '''
    load meshmask, expand fields
    just those needed to degrade forcings
    '''
    M = xr.open_dataset(maskfile)
    M1 = {}
    #
    M1["e1t"] = dm.xpnd_wrap(M["e1t"], 'interp', ndeg)
    M1["e1u"] = dm.xpnd_wrap(M["e1u"], 'interp', ndeg)
    M1["e1v"] = dm.xpnd_wrap(M["e1v"], 'interp', ndeg)
    M1["e2t"] = dm.xpnd_wrap(M["e2t"], 'edge', ndeg)
    M1["e2u"] = dm.xpnd_wrap(M["e2u"], 'edge', ndeg)
    M1["e2v"] = dm.xpnd_wrap(M["e2v"], 'edge', ndeg)
    M1["e3t_0"] = dm.xpnd_wrap(M["e3t_0"], 'edge', ndeg)
    M1["e3u_0"] = dm.xpnd_wrap(M["e3u_0"], 'edge', ndeg)
    M1["e3v_0"] = dm.xpnd_wrap(M["e3v_0"], 'edge', ndeg)
    M1["glamt"] = dm.xpnd_wrap(M["glamt"], 'interp', ndeg)
    M1["glamu"] = dm.xpnd_wrap(M["glamu"], 'interp', ndeg)
    M1["glamv"] = dm.xpnd_wrap(M["glamv"], 'interp', ndeg)
    M1["gphit"] = dm.xpnd_wrap(M["gphit"], 'interp', ndeg)
    M1["gphiu"] = dm.xpnd_wrap(M["gphiu"], 'interp', ndeg)
    M1["gphiv"] = dm.xpnd_wrap(M["gphiv"], 'interp', ndeg)
    M1["tmask"] = dm.xpnd_wrap(M["tmask"], 'edge', ndeg) #there is some error here!
    M1["umask"] = dm.xpnd_wrap(M["umask"], 'edge', ndeg) #there is some error here!
    M1["vmask"] = dm.xpnd_wrap(M["vmask"], 'edge', ndeg) #there is some error here!
    # Bonus variable for degrading forcings
    M1['h_column_t'] = M1['e3t_0'].sum(dim='z')
    #
    M.close()
    M1 = xr.Dataset(M1)
    return M1


def load_coords_degraded(maskfile_d):
    '''
    loads lat, lon coordinates on T, U, V grids
    from already degraded mask
    to avoid having to degrade each time the same arrays
    need to run degrade_mesh.py first!
    '''
    Md = xr.open_dataset(maskfile_d)
    C = Md[["glamt", "glamu", "glamv", "gphit", "gphiu", "gphiv"]]
    C['glamw'] = C['glamt']
    C['gphiw'] = C['gphit']
    return C

def load_tfile(infile, ndeg=1):
    '''
    loads t-grid forcing file, pads a rim all around it
    to make j, i multiples of ndeg,
    populates rim with edge values
    (ugly, but it's the atlantic buffer that gets populated,
    the model won't use those values anyway)
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    # most of these ogstm does not need (commented)
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['votemper'] = dm.xpnd_wrap(F['votemper'], 'edge', ndeg)
    F1['vosaline'] = dm.xpnd_wrap(F['vosaline'], 'edge', ndeg)
    F1['sossheig'] = dm.xpnd_wrap(F['sossheig'], 'edge', ndeg)
    #F1['sossh_ib'] = dm.xpnd_wrap(F['sossh_ib'], 'edge', ndeg)
    #F1['sowaflup'] = dm.xpnd_wrap(F['sowaflup'], 'edge', ndeg)
    #F1['soevapor'] = dm.xpnd_wrap(F['soevapor'], 'edge', ndeg)
    #F1['soprecip'] = dm.xpnd_wrap(F['soprecip'], 'edge', ndeg)
    F1['sorunoff'] = dm.xpnd_wrap(F['sorunoff'], 'edge', ndeg)
    F1['soshfldo'] = dm.xpnd_wrap(F['soshfldo'], 'edge', ndeg)
    #F1['sohefldo'] = dm.xpnd_wrap(F['sohefldo'], 'edge', ndeg)
    #F1['solofldo'] = dm.xpnd_wrap(F['solofldo'], 'edge', ndeg)
    #F1['sosefldo'] = dm.xpnd_wrap(F['sosefldo'], 'edge', ndeg)
    #F1['solafldo'] = dm.xpnd_wrap(F['solafldo'], 'edge', ndeg)
    F1['somxl010'] = dm.xpnd_wrap(F['somxl010'], 'edge', ndeg)
    # overwrite coordinates, else it complains
    for vv in F1.keys():
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    #
    F1['deptht'] = F['deptht']
    F1['deptht_bounds'] = F['deptht_bounds']
    F1['time_instant'] = F['time_instant']
    F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_ufile(infile, ndeg=1):
    '''
    as in load_tfile()
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vozocrtx'] = dm.xpnd_wrap(F['vozocrtx'], 'edge', ndeg)
    F1['sozotaux'] = dm.xpnd_wrap(F['sozotaux'], 'edge', ndeg)
    vlist = list(F1.keys())
    # overwrite coordinates, else it complains
    for vv in vlist:
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthu'] = F['depthu'] 
    F1['depthu_bounds'] = F['depthu_bounds'] 
    F1['time_instant'] = F['time_instant'] 
    F1['time_instant_bounds'] = F['time_instant_bounds'] 
    F1['time_counter'] = F['time_counter'] 
    F1['time_counter_bounds'] = F['time_counter_bounds'] 
    F1['time_centered'] = F['time_centered'] 
    F1['time_centered_bounds'] = F['time_centered_bounds'] 
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_vfile(infile, ndeg=1):
    '''
    as in load_tfile(), load_ufile
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vomecrty'] = dm.xpnd_wrap(F['vomecrty'], 'edge', ndeg)
    F1['sometauy'] = dm.xpnd_wrap(F['sometauy'], 'edge', ndeg)
    vlist = list(F1.keys())
    # overwrite coordinates, else it complains
    for vv in vlist:
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthv'] = F['depthv']
    F1['depthv_bounds'] = F['depthv_bounds']
    F1['time_instant'] = F['time_instant']
    F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_wfile(infile, ndeg=1):
    '''
    as in load_tfile(), load_ufile(), load vfile()
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vovecrtz'] = dm.xpnd_wrap(F['vovecrtz'], 'edge', ndeg)
    F1['votkeavt'] = dm.xpnd_wrap(F['votkeavt'], 'edge', ndeg)
    # overwrite coordinates, else it complains
    for vv in F1.keys():
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthw'] = F['depthw']
    F1['depthw_bounds'] = F['depthw_bounds']
    #F1['time_instant'] = F['time_instant']
    #F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1


def load_fffile(infile, tuv, ndeg=1):
    if tuv == 'T':
        F1 = load_tfile(infile, ndeg)
    elif tuv == 'U':
        F1 = load_ufile(infile, ndeg)
    elif tuv == 'V':
        F1 = load_vfile(infile, ndeg)
    elif tuv == 'W':
        F1 = load_wfile(infile, ndeg)
    return F1

def get_mask_fields(M):
    '''
    gets the arrays needed to compute cell surface areas and volumes
    those from the mask that do not change
    so it doesn't have to retrieve them each time
    '''
    h_column_t = M['h_column_t'].values[:]
    tmask =           M['tmask'].values[:]
    umask =           M['umask'].values[:]
    vmask =           M['vmask'].values[:]
    e3t_0 =           M['e3t_0'].values[:]
    e3u_0 =           M['e3u_0'].values[:]
    e3v_0 =           M['e3v_0'].values[:]
    e1u =               M['e1u'].values[:]
    e2u =               M['e2u'].values[:]
    e1v =               M['e1v'].values[:]
    e2v =               M['e2v'].values[:]
    e1t =               M['e1t'].values[:]
    e2t =               M['e2t'].values[:]
    At = e1t * e2t #e1v * e2u
    Aw = np.repeat(At, repeats=tmask.shape[1], axis=1)
    return h_column_t, tmask, umask, vmask, e3t_0, e3u_0, e3v_0, e1u, e2u, e1v, e2v, e1t, e2t, At, Aw

def get_nanmask(mask):
    nanmask = np.ones(mask.shape)
    nanmask[mask==0.0] = np.nan
    return nanmask

def get_weights(M, T):
    '''
    given the meshmask and the current t-grid file (with free surface)
    computes cell surface areas on U, V grids and cell volume on T grid
    and surface area on T / W grid
    algorithm from ogstm/src/IO/forcing_phys.f90
    (NB, this one pure numpy because xarray is not good at broadcasting)
    '''
    # THIS GOES OUTSIDE THE FUNCTION, SO I DON'T DO IT EACH TIME I CALL IT
    h_column_t, tmask, umask, vmask, e3t_0, e3u_0, e3v_0, e1u, e2u, e1v, e2v, e1t, e2t, At, Aw = get_mask_fields(M)
    nan_tmask = get_nanmask(tmask) 
    nan_umask = get_nanmask(umask)
    nan_vmask = get_nanmask(vmask)
    At = At * nan_tmask[0,0,:,:] #used for degrading T-grid, needs NaN (weight only valid values)
    # THIS GOES OUTSIDE THE FUNCTION, SO I DON'T DO IT EACH TIME I CALL IT
    ssh = tmask * T['sossheig'].values[:]
    correction_e3t = 1.0 + (ssh / h_column_t)
    e3t = tmask * e3t_0 * correction_e3t
    #
    e1u_x_e2u = (e1u * e2u)
    e1v_x_e2v = (e1v * e2v)
    e1t_x_e2t = (e1t * e2t)
    diff_e3t = e3t - e3t_0
    #
    s0 = e1t_x_e2t * diff_e3t
    s1 = np.zeros(s0.shape)
    s2 = np.zeros(s0.shape)
    s1[:,:,:,1:] = e1t_x_e2t[:,:,:,1:] * diff_e3t[:,:,:,1:] #THIS INDEXING SHOULD BE CORRECT
    s2[:,:,1:,:] = e1t_x_e2t[:,:,1:,:] * diff_e3t[:,:,1:,:] #THIS INDEXING SHOULD BE CORRECT
    #
    e3u = e3u_0 + (0.5 * (umask / e1u_x_e2u) * (s0 + s1))
    e3v = e3v_0 + (0.5 * (vmask / e1v_x_e2v) * (s0 + s2))
    #
    B = {}
    B['e1v'] = xr.DataArray(e1v, name='e1v' , dims=('time','z_a','y','x'))
    B['e2u'] = xr.DataArray(e2u, name='e2u' , dims=('time','z_a','y','x'))
    B['At'] = xr.DataArray(At, name='At' , dims=('time','z_a','y','x'))
    B['Aw'] = xr.DataArray(Aw, name='Aw' , dims=('time','z','y','x'))
    B['Au'] = xr.DataArray((e2u * e3u), name='Au' , dims=('time','z','y','x'))
    B['Av'] = xr.DataArray((e1v * e3v), name='Av' , dims=('time','z','y','x'))
    B['V'] = xr.DataArray((nan_tmask * e1v * e2u * e3t), name='V' , dims=('time','z','y','x')) #for degrading T-grid, don't weight land values
    B = xr.Dataset(B)
    return B

def degrade_wrap(X, degr_op, B=1.0, ndeg=1):
    print('---1---')
    X = dm.reshape_blocks(X, ndeg)
    print('---2---')
    X = degr_op(X, ndeg, B)
    print('...')
    return X

def init_ds(D, C, tuv):
    '''
    adds to degraded dataset the variables that are always the same:
    nav_lat, nav_lon, depth, 
    and those that do not need regridding, time
    '''
    Dd = {}
    Dd[f'nav_lon'] = C[f'glam{tuv}']
    Dd[f'nav_lat'] = C[f'gphi{tuv}']
    Dd[f'depth{tuv}'] = D[f'depth{tuv}']
    Dd[f'depth{tuv}_bounds'] = D[f'depth{tuv}_bounds']
    #Dd['time_instant'] = D['time_instant']
    Dd['time_counter'] = D['time_counter']
    Dd['time_centered'] = D['time_centered']
    #Dd['time_instant_bounds'] = D['time_instant_bounds']
    Dd['time_counter_bounds'] = D['time_counter_bounds']
    Dd['time_centered_bounds'] = D['time_centered_bounds']
    if tuv != 'w':
        # for some reason these are not in W-grid files
        Dd['time_instant'] = D['time_instant']
        Dd['time_instant_bounds'] = D['time_instant_bounds']
    return Dd

def degrade_V(V, B, C, ndeg=1):
    '''
    degrades resolution of V-grid file
    vomecrty
    sometauy
    '''
    B = B.rename({'time':'time_counter', 'z':'depthv'})
    #Vd = {}
    Vd = init_ds(V, C, 'v')
    #
    Vd['vomecrty'] = degrade_wrap(V['vomecrty'], dm.vawmean_jstep, B['Av'], ndeg)
    Vd['sometauy'] = degrade_wrap(V['sometauy'], dm.e1vwmean_jstep, B['e1v'], ndeg)
    #
    Vd = xr.Dataset(Vd)
    return Vd

def degrade_U(U, B, C, ndeg=1):
    '''
    degrades resolution of U-grid file
    vozocrtx
    sozotaux
    '''
    B = B.rename({'time':'time_counter', 'z':'depthu'})
    #Ud = {}
    Ud = init_ds(U, C, 'u')
    #
    Ud['vozocrtx'] = degrade_wrap(U['vozocrtx'], dm.uawmean_istep, B['Au'], ndeg)
    Ud['sozotaux'] = degrade_wrap(U['sozotaux'], dm.e2uwmean_istep, B['e2u'], ndeg)
    #
    Ud = xr.Dataset(Ud)
    return Ud

def degrade_W(W, B, C, ndeg=1):
    '''
    degrades resolution of W-grid file
    vovecrtz
    votkeavt
    '''
    B = B.rename({'time':'time_counter', 'z':'depthw'})
    #Wd = {}
    Wd = init_ds(W, C, 'w')
    #
    # use dm.vwmean(), but actually uses 3D Aw as weight
    Wd['vovecrtz'] = degrade_wrap(W['vovecrtz'], dm.vwmean, B['Aw'], ndeg)
    Wd['votkeavt'] = degrade_wrap(W['votkeavt'], dm.vwmean, B['Aw'], ndeg)
    #
    Wd = xr.Dataset(Wd)
    return Wd


def degrade_T(T, B, C, ndeg=1):
    '''
    degrades resolution of T-grid file
    '''
    B = B.rename({'time':'time_counter', 'z':'deptht'})
    Td = init_ds(T, C, 't')
    #Td = {}
    #
    Td['votemper'] = degrade_wrap(T['votemper'], dm.vwmean, B['V'], ndeg) 
    Td['vosaline'] = degrade_wrap(T['vosaline'], dm.vwmean, B['V'], ndeg) 
    Td['sossheig'] = degrade_wrap(T['sossheig'], dm.awmean, B['At'], ndeg) 
    Td['soshfldo'] = degrade_wrap(T['soshfldo'], dm.awmean, B['At'], ndeg) 
    Td['sorunoff'] = degrade_wrap(T['sorunoff'], dm.awmean, B['At'], ndeg)
    Td['somxl010'] = degrade_wrap(T['somxl010'], dm.awmean, B['At'], ndeg)
    #
    Td = xr.Dataset(Td)
    return Td

def make_outdir(outdir, outfile):
    #/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/qDEG_SETUP/FORCINGS/
    #T20020308-12:00:00.nc
    yyyy = outfile[1:5]
    mm = outfile[5:7]
    outdir = f'{outdir}/{yyyy}/{mm}/'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True) 
    else:
        pass
    return outdir

def degrade(F, tuv, B, C, outdir, outfile, ndeg=1):
    if tuv == 'T':
        Fd = degrade_T(F, B, C, ndeg)
    elif tuv == 'U':
        Fd = degrade_U(F, B, C, ndeg)
    elif tuv == 'V':
        Fd = degrade_V(F, B, C, ndeg)
    elif tuv == 'W':
        Fd = degrade_W(F, B, C, ndeg)
    outdir = make_outdir(outdir, outfile)
    outfile = outdir+outfile
    dump_netcdf(Fd, outfile)
    return Fd

def get_flist(tuvw, Params):
    ffdir = Params['ffdir']
    y0 = Params['y0']
    yE = Params['yE']
    flist = [glob(f'{ffdir}/{YYYY}/??/{tuvw}*.nc') for YYYY in range(y0, yE+1)]
    flist = sorted(list(chain.from_iterable(flist)))
    return flist

if __name__=='__main__':
    #
    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nranks = comm.size
    except:
        comm = None
        rank = 0
        nranks = 1
    print (rank,nranks)
    #
    Params = load_parameters()
    #
    maskfile = Params['maskfile']
    maskfile_d = Params['maskfile_d']
    ffdir = Params['ffdir']
    outdir = Params['outdir']
    y0 = Params['y0']
    yE = Params['yE']
    #
    ndeg = Params['ndeg']
    #
    flistt = get_flist('T', Params) 
    flistu = get_flist('U', Params)
    flistv = get_flist('V', Params)
    flistw = get_flist('W', Params)
    #
    # load mesh and expand lat, lon to make them multiples of ndeg
    if rank==0:
        print('loading mask')
        M = load_mesh_light(maskfile, ndeg) #NB, here Au, Av, V, are calculated with e3tuv_0
        C = load_coords_degraded(maskfile_d)
    else:
        M = None
        C = None
    #
    if nranks > 1:
        M = comm.bcast(M, root=0)
        C = comm.bcast(C, root=0)
    #
    # LOOP ON FILES
    print('loading T, U, V files')
    ziter = list(zip(flistt, flistu, flistv, flistw))
    for tfile, ufile, vfile, wfile in ziter[rank::nranks]:
        V = load_fffile(vfile, 'V', ndeg)
        U = load_fffile(ufile, 'U', ndeg)
        W = load_fffile(wfile, 'W', ndeg)
        T = load_fffile(tfile, 'T', ndeg)
        B = get_weights(M, T)
        # degrade mesh
        print('degrado')
        outft = tfile.split('/')[-1]
        outfu = ufile.split('/')[-1]
        outfv = vfile.split('/')[-1]
        outfw = wfile.split('/')[-1]
        degrade(T, 'T', B, C, outdir, outft, ndeg)
        degrade(U, 'U', B, C, outdir, outfu, ndeg)
        degrade(V, 'V', B, C, outdir, outfv, ndeg)
        degrade(W, 'W', B, C, outdir, outfw, ndeg)


