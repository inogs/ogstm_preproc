import numpy as np
import xarray as xr
from argparse import ArgumentParser
import yaml

'''
generates a reduced horizontal resolution mesh from an original mesh (with Atlantic buffer)
vertical resolution is unchanged.
resolution is degraded by merging cells in (ndeg x ndeg) blocks.
functining: 
0. load original mesh, 
1. expand lat, lon to make the new size a multiple of ndeg, extrapolate values properly
2. compute cell surface areas and volume (for weighted means)
3. degrade resolution of original variables
4. dump new meshmask.nc
usage: python degrade_mesh.py -y degrade_mesh.yaml

yaml file example:
maskfile: '/path/to/meshmask_in.nc'
outfile: '/path/to/meshmask_out.nc'
outfile_med: '/path/to/meshmask_out_med.nc'
ndeg: 6    #nr of cells to join in i, j
thresh: 1  #water point per block to assign to waterpoint in destination grid

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
    # yamlfile = 'degrade_mesh.yaml'
    with open(yamlfile) as f:
        Params = yaml.load(f, Loader=yaml.Loader)
    return(Params)


# LOAD MESHMASK

def adjust_dims(X, n=4):
    '''
    add singleton dimensions in order to have it 4D
    (t, z, y, x)
    '''
    cdim = 0
    nd = X.ndim
    while X.ndim < n:
        new_dim = f'j{cdim}'
        X = X.expand_dims({new_dim: 1})
        cdim = cdim + 1
    return X, nd
#
def deadjust_dims(X, nd):
    '''
    remove some singleton dimensions 
    in order to get back to original dimensions
    '''
    while X.ndim > nd and X.shape[0] == 1:
        X = X.squeeze(axis=0)
    return X

def fill_rim_lin_2D(X):
    '''
    fills an outer rim of NaN by linear extrapolation of internal values
    in x and y direction (works also with constant values)
    for glamt, gphit, etc.
    '''
    t, w, j, i = X.shape
    # detect rim
    X2d = X.squeeze().values
    isnan = np.isnan(X2d)
    south_rim = np.argmax((~isnan).any(axis=1))       # first valid row
    north_rim = np.argmax((~isnan[::-1]).any(axis=1)) # last valid row from top
    north_rim = j - north_rim                         # last valid row
    west_rim  = 0                                     # first valid column (Atl open boundary)
    east_rim  = np.argmin((~isnan[::-1]).any(axis=0)) # last valid colum
    # linear interpolation
    filled = X2d.copy()
    dy = filled[north_rim-1,:] - filled[north_rim-2,:]
    # fill south
    filled[:south_rim, :] = filled[south_rim,:] - np.outer(np.arange(south_rim)[::-1]+1, dy)
    # fill north
    filled[north_rim:, :] = filled[north_rim-1,:] + np.outer(np.arange(j-north_rim)+1, dy)
    # fill west
    dx = filled[:,1] - filled[:,0]
    filled[:, :west_rim] = filled[:, west_rim][:,None] - np.outer(dx, np.arange(west_rim)[::-1]+1)
    # fill east
    filled[:, east_rim:] = filled[:, east_rim-1][:,None] + np.outer(dx, np.arange(i-east_rim))
    # out
    X[...,:,:] = filled
    return X

def expand_array(X, pad_op, ndeg=1):
    '''
    pad array with linear extrapolation (pad_op='interp')
    or edge vaslues (pad_op='edge')
    to make j, i multiples of ndeg.
    also add a rim of size ndeg, so I'm sure to have
    land points at the N, S, E boundaries
    when degrading resolution
    '''
    t, w, j, i = X.shape
    dj = (ndeg - j % ndeg) % ndeg
    di = (ndeg - i % ndeg) % ndeg
    pw = {'y':(dj+ndeg, ndeg), 'x':(0, di+ndeg)}
    if pad_op=='edge':
        X = X.pad(pad_width=pw, mode='edge')
    elif pad_op=='interp':
        X = X.pad(pad_width=pw, mode='constant', constant_values=np.nan)
        X = fill_rim_lin_2D(X)
    return X

def xpnd_wrap(X, pad_op, ndeg=1):
    X, nd = adjust_dims(X)
    X = expand_array(X, pad_op, ndeg)
    X = deadjust_dims(X, nd)
    return X

def e3_2_gdep(M, tuvw):
    '''
    reconstructs gdept from e3t_0
    cause no gdept in the original Nemo mesh
    '''
    e3 = M[f'e3{tuvw}_0']
    gdep = e3.cumsum(dim="z") - (e3 / 2.0)
    return gdep

def load_mesh(maskfile, ndeg=1):
    '''
    load meshmask, expand fields
    '''
    M = xr.open_dataset(maskfile)
    M1 = {}
    #
    M1["coastp"] = xpnd_wrap(M["coastp"], 'edge', ndeg)
    M1["e1f"] = xpnd_wrap(M["e1f"], 'interp', ndeg)
    M1["e1t"] = xpnd_wrap(M["e1t"], 'interp', ndeg)
    M1["e1u"] = xpnd_wrap(M["e1u"], 'interp', ndeg)
    M1["e1v"] = xpnd_wrap(M["e1v"], 'interp', ndeg)
    M1["e2f"] = xpnd_wrap(M["e2f"], 'edge', ndeg)
    M1["e2t"] = xpnd_wrap(M["e2t"], 'edge', ndeg)
    M1["e2u"] = xpnd_wrap(M["e2u"], 'edge', ndeg)
    M1["e2v"] = xpnd_wrap(M["e2v"], 'edge', ndeg)
    M1["e3t_0"] = xpnd_wrap(M["e3t_0"], 'edge', ndeg)
    M1["e3u_0"] = xpnd_wrap(M["e3u_0"], 'edge', ndeg)
    M1["e3v_0"] = xpnd_wrap(M["e3v_0"], 'edge', ndeg)
    M1["e3w_0"] = xpnd_wrap(M["e3w_0"], 'edge', ndeg)
    M1["ff"] = xpnd_wrap(M["ff"], 'interp', ndeg)
    M1["gdept"] = xpnd_wrap(e3_2_gdep(M, 't'), 'edge', ndeg) #M["gdept"]
    M1["gdepw"] = xpnd_wrap(e3_2_gdep(M, 'w'), 'edge', ndeg) #M["gdepw"]
    M1["glamf"] = xpnd_wrap(M["glamf"], 'interp', ndeg)
    M1["glamt"] = xpnd_wrap(M["glamt"], 'interp', ndeg)
    M1["glamu"] = xpnd_wrap(M["glamu"], 'interp', ndeg)
    M1["glamv"] = xpnd_wrap(M["glamv"], 'interp', ndeg)
    M1["gphif"] = xpnd_wrap(M["gphif"], 'interp', ndeg)
    M1["gphit"] = xpnd_wrap(M["gphit"], 'interp', ndeg)
    M1["gphiu"] = xpnd_wrap(M["gphiu"], 'interp', ndeg)
    M1["gphiv"] = xpnd_wrap(M["gphiv"], 'interp', ndeg)
    M1["nav_lat"] = xpnd_wrap(M["nav_lat"], 'interp', ndeg)
    M1["nav_lon"] = xpnd_wrap(M["nav_lon"], 'interp', ndeg)
    M1["nav_lev"] = M["nav_lev"]
    M1["fmask"] = xpnd_wrap(M["fmask"], 'edge', ndeg) #there is some error here!
    M1["tmask"] = xpnd_wrap(M["tmask"], 'edge', ndeg) #there is some error here!
    M1["umask"] = xpnd_wrap(M["umask"], 'edge', ndeg) #there is some error here!
    M1["vmask"] = xpnd_wrap(M["vmask"], 'edge', ndeg) #there is some error here!
    # Compute cell surface areas and volume
    M1['At'] = xr.DataArray(data = M1['e1t'].values[:] * M1['e2t'].values[:], 
                            dims = ('time','z_a','y','x'),
                            name = 'At')
    M1['Au'] = xr.DataArray(data = M1['e2u'].values[:] * M1['e3u_0'].values[:],
                            dims = ('time','z','y','x'),
                            name = 'Au')
    M1['Av'] = xr.DataArray(data = M1['e2v'].values[:] * M1['e3v_0'].values[:],
                            dims = ('time','z','y','x'),
                            name = 'Av')
    M1['V'] = xr.DataArray(data = M1['e1v'].values[:] * M1['e2u'].values[:] * M1['e3t_0'].values[:],
                            dims = ('time','z','y','x'),
                            name = 'V')
    # Bonus variable for degrading forcings
    M1['h_column_t'] = (M1['tmask'] * M1['e3t_0']).sum(dim='z')
    #
    M.close()
    M1 = xr.Dataset(M1)
    return M1

# DEGRADE RESOLUTION

def reshape_blocks(X, ndeg=1):
    '''
    reshape in blocks of shape ndeg along j, i
    xarray's convoluted way of doing np.reshape,
    assume I'm always reshaping dims called y, x
    chunking is necessary to avoid OOM
    (still cuz xarray's way of reshaping is not super smart)
    '''
    cs = ndeg * 5
    X = X.chunk({'y':cs, 'x':cs})
    X = X.coarsen(y=ndeg, x=ndeg).construct(y=('y', 'y_b'), x=('x', 'x_b'))

    # X has the dimension y_b, x_b of size ndeg
    return X

# REGRIDDING FUNCTIONS
def isum_jmean(X, ndeg=1, W=1.0):
    '''
    sum along longitude (i) & mean along latitude (j)
    no need to weight by area! for e1t
    '''
    X = X.sum(dim='x_b', skipna=True)  #blocks sum along lon (i)
    X = X.mean(dim='y_b', skipna=True) #block mean along lat (j)
    return X

def imean_jstep(X, ndeg=1, W=1.0):
    '''
    mean along i (lon)
    & one element each ndeg along j (lat)
    glamt
    '''
    X = X.isel({'y_b':-1})
    X = X.mean(dim='x_b', skipna=True)
    return X

def jmean_istep(X, ndeg=1, W=1.0):
    '''
    mean along j (lat)
    & one element each ndeg along i (lon)
    gphit
    '''
    X = X.isel({'x_b':-1})
    X = X.mean(dim='y_b', skipna=True)
    return X

def isum_jmean(X, ndeg=1, W=1.0):
    '''
    sum along longitude (i) & mean along latitude (j)
    no need to weight by area!
    e1t
    '''
    X = X.sum(dim='x_b', skipna=True)  #blocks sum along lon (i)
    X = X.mean(dim='y_b', skipna=True) #block mean along lat (j)
    return X

def jsum_imean(X, ndeg=1, W=1.0):
    '''
    sum along latitude (j) & mean along longitude (i)
    no need to weight by area!
    e2t
    '''
    X = X.sum(dim='y_b', skipna=True)  #blocks sum along lat (j)
    X = X.mean(dim='x_b', skipna=True) #block mean along lon (i)
    return X

def isum_jstep(X, ndeg=1, W=1.0):
    '''
    sum along longitude (i) & one point each ndeg along latitude (j)
    '''
    X = X.sum(dim='x_b', skipna=True)
    X = X.isel({'y_b': -1})
    return X

def jsum_istep(X, ndeg=1, W=1.0):
    '''
    sum along latitude (j) & one point each ndeg along longitude (i)
    '''
    X = X.sum(dim='y_b', skipna=True)
    X = X.isel({'x_b': -1})
    return X

def awmean(X, ndeg=1, W=1.0):
    '''
    mean, weighted by surface area (t-grid)
    '''
    At = reshape_blocks(W, ndeg)
    X = (At * X).sum(dim=('y_b', 'x_b'), skipna=True) / At.sum(dim=('y_b', 'x_b'), skipna=True)
    return X

def vwmean(X, ndeg=1, W=1.0):
    '''
    W is a dataarray with cell volumes 3D
    mean, weighted by surface area (t-grid)
    '''
    V = reshape_blocks(W, ndeg)
    X = (V * X).sum(dim=('y_b', 'x_b'), skipna=True) / V.sum(dim=('y_b', 'x_b'), skipna=True)
    return X

def uawmean_istep(X, ndeg=1, W=1.0):
    '''
    mean, weighted by lateral u-grid area along lat (j)
    & one element each ndeg along lon (i)
    '''
    Au = reshape_blocks(W, ndeg)
    Au = Au.isel({'x_b':-1})
    X = X.isel({'x_b':-1})
    X = (Au * X).sum(dim='y_b', skipna=True) / Au.sum(dim='y_b', skipna=True)
    return X

def vawmean_jstep(X, ndeg=1, W=1.0):
    '''
    mean, weighted by lateral u-grid area along lat (j)
    & one element each ndeg along lon (i)
    '''
    Av = reshape_blocks(W, ndeg)
    Av = Av.isel({'y_b':-1})
    X = X.isel({'y_b':-1})
    X = (Av * X).sum(dim='x_b', skipna=True) / Av.sum(dim='x_b', skipna=True)
    return X

def e1vwmean_jstep(X, ndeg=1, W=1.0):
    '''
    mean, weighted by e1v, v-grid length along lon (i)
    & one element each ndeg along lat (j)
    sometauy (from phy forcings)
    '''
    e1v = reshape_blocks(W, ndeg)
    e1v = e1v.isel({'y_b':-1})
    X = X.isel({'y_b':-1})
    X = (e1v * X).sum(dim='x_b', skipna=True) / e1v.sum(dim='x_b', skipna=True)
    return X

def e2uwmean_istep(X, ndeg=1, W=1.0):
    '''
    mean, weighted by e1u, u-grid length along lat (j)
    & one element each ndeg along lon (i)
    sozotaux (from phy forcings)
    '''
    e2u = reshape_blocks(W, ndeg)
    e2u = e2u.isel({'x_b':-1})
    X = X.isel({'x_b':-1})
    X = (e2u * X).sum(dim='y_b', skipna=True) / e2u.sum(dim='y_b', skipna=True)
    return X

def waterpt_thresh(X24, thresh=1 ,ndeg=1):
    '''
    Args:
    X: tmask
    1 or 0 if enough water points (or volume?)
    the higher thresh, the least water points needed 
    to classify a cell as water
    thresh=ndeg**2: all 36 points must be water
    thresh=1 (36) just one point needed
    tmask
    '''
    X24 = reshape_blocks(X24, ndeg)
    thr = ndeg**2 - thresh + 1
    X_degr = X24.sum(dim=('x_b', 'y_b'))
    X_degr= (X_degr >= (ndeg**2 // thr)).astype(int)
    return X_degr

def from_tmask(tmask):
    '''
    computes umask, vmask, fmask from tmask
    after Madec 2015 NEMO Manual, pp. 67
    NB Madec's fmask definition is not correct!
    tmask.shape: (1, 141, 380, 1307)
    '''
    # CHECK I'M PADDING AT THE RIGHT END OF THE ARRAY!!!
    tmask_ipad = tmask.shift(x=1, fill_value=0.0)
    tmask_jpad = tmask.shift(y=1, fill_value=0.0)
    tmask_ijpad = tmask.shift(x=1, y=1, fill_value=0.0)
    #
    umask = tmask * tmask_ipad
    vmask = tmask * tmask_jpad
    # ERROR: THIS fmask DOES NOT CORRESPOND WITH WHAT'S IN meshmask.nc
    fmask = tmask * tmask_ipad * tmask_jpad * tmask_ijpad
    #
    return umask, vmask, fmask

def noop(X):
    '''
    no operation, for variables that need no regridding
    e.g. nav_lev
    '''
    return X

# END REGRIDDING OPERATIONS

def degr_wrap(X24, degr_op, ndeg=1, W=1.0):
    X24, nd = adjust_dims(X24) # make 4D
    X24_dims = reshape_blocks(X24, ndeg)
    X_degr = degr_op(X24_dims, ndeg, W)
    X_degr = deadjust_dims(X_degr, nd) # back to original dims
    return X_degr

def degrade_mesh(M1, thresh=1, ndeg=1):
    '''
    generates new reduced mesh for each variable
    '''
    M2 = {}
    #
    M2["coastp"] = degr_wrap(M1["coastp"], awmean, ndeg, W=M1['At'])
    M2["e1f"] = degr_wrap(M1["e1f"], isum_jstep, ndeg, W=None)
    M2["e1t"] = degr_wrap(M1["e1t"], isum_jmean, ndeg, W=None)
    M2["e1u"] = degr_wrap(M1["e1u"], isum_jmean, ndeg, W=None)
    M2["e1v"] = degr_wrap(M1["e1v"], isum_jstep, ndeg, W=None)
    M2["e2f"] = degr_wrap(M1["e2f"], jsum_istep, ndeg, W=None)
    M2["e2t"] = degr_wrap(M1["e2t"], jsum_istep, ndeg, W=None)
    M2["e2u"] = degr_wrap(M1["e2u"], jsum_istep, ndeg, W=None)
    M2["e2v"] = degr_wrap(M1["e2v"], jsum_istep, ndeg, W=None)
    M2["e3t_0"] = degr_wrap(M1["e3t_0"], vwmean, ndeg, W=M1['V'])
    M2["e3u_0"] = degr_wrap(M1["e3u_0"], uawmean_istep, ndeg, W=M1['Au'])
    M2["e3v_0"] = degr_wrap(M1["e3v_0"], vawmean_jstep, ndeg, W=M1['Av'])
    M2["e3w_0"] = degr_wrap(M1["e3w_0"], vwmean, ndeg, W=M1['V'])
    M2["ff"] = degr_wrap(M1["ff"], jmean_istep, ndeg, W=None)
    #M2["gdept"] = noop(M1["gdept"])
    #M2["gdepw"] = noop(M1["gdepw"])
    M2["gdept"] = degr_wrap(M1["gdept"], vwmean, ndeg, W=M1['V'])
    M2["gdepw"] = degr_wrap(M1["gdepw"], vwmean, ndeg, W=M1['V'])
    M2["glamf"] = degr_wrap(M1["glamf"], jmean_istep, ndeg, W=None)
    M2["glamt"] = degr_wrap(M1["glamt"], imean_jstep, ndeg, W=None)
    M2["glamu"] = degr_wrap(M1["glamu"], jmean_istep, ndeg, W=None)
    M2["glamv"] = degr_wrap(M1["glamv"], imean_jstep, ndeg, W=None)
    M2["gphif"] = degr_wrap(M1["gphif"], imean_jstep, ndeg, W=None)
    M2["gphit"] = degr_wrap(M1["gphit"], jmean_istep, ndeg, W=None)
    M2["gphiu"] = degr_wrap(M1["gphiu"], jmean_istep, ndeg, W=None)
    M2["gphiv"] = degr_wrap(M1["gphiv"], imean_jstep, ndeg, W=None)
    M2["nav_lat"] = degr_wrap(M1["nav_lat"], jmean_istep, ndeg, W=None)
    M2["nav_lon"] = degr_wrap(M1["nav_lon"], imean_jstep, ndeg, W=None)
    M2["nav_lev"] = noop(M1["nav_lev"])
    M2["tmask"] = waterpt_thresh(M1["tmask"], thresh ,ndeg)
    umask, vmask, fmask = from_tmask(M2["tmask"])
    M2["umask"] = umask
    M2["vmask"] = vmask
    M2["fmask"] = fmask
    #
    M2 = xr.Dataset(M2)
    return M2

def cut_med(M2, lon_cut=-8.875, depth_cut=4153.0, biscay_land=True):
    '''
    cut out atlantic buffer,
    mask bay of biscay
    '''
    lon1D = M2['glamt'].isel(time=0, z_a=0, y=0)
    z1D = M2['nav_lev']
    #
    x_sel = np.arange(len(lon1D))[lon1D >= lon_cut]
    z_sel = np.arange(len(z1D))[z1D <= depth_cut]
    #
    MMed = M2.isel(z=z_sel, x=x_sel)
    if biscay_land:
        bool1 = (MMed.glamt < 0) & (MMed.gphit > 42)             
        bool2 = (MMed.glamt < -6.) & (MMed.gphit > 37.25)
        land = (bool1 | bool2).squeeze()
        MMed['tmask'] = MMed['tmask'] * ~land
        MMed['umask'] = MMed['umask'] * ~land
        MMed['vmask'] = MMed['vmask'] * ~land
        MMed['fmask'] = MMed['fmask'] * ~land
    return MMed

def dump_mesh(M2, outfile):
    print(f'saving: {outfile}')
    M2.to_netcdf(outfile, mode='w', unlimited_dims='time', format='NETCDF4_CLASSIC')
    return

if __name__=='__main__':
    Params = load_parameters()
    maskfile = Params['maskfile']
    outfile = Params['outfile']
    outfile_med = Params['outfile_med']
    # by how much cells to degrade (nr of cells to join in i, j)
    ndeg = Params['ndeg']
    # threshold of points to assign waterpoint to degraded mesh
    thresh = Params['thresh'] #THE ONLY OPTION THAT ENSURES NO RIVERS END UP ON LAND
    # load meshmask and expand dimensions
    print('---1---')
    M1 = load_mesh(maskfile, ndeg)
    #
    print('---2---')
    M2 = degrade_mesh(M1, thresh, ndeg)
    #
    print('---3---')
    dump_mesh(M2, outfile)
    #
    print('---4---')
    MMed = cut_med(M2)
    dump_mesh(MMed, outfile_med)


