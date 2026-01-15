import numpy as np
import xarray as xr
import regridding as rg
from commons import reshape_blocks

def adjust_dims(X: xr.DataArray, n: int = 4):
    '''
    add singleton dimensions in front of existing ones in order to have n dimensions
    
    Args:
    X: xarray DataArray, with dims < = n
    n: desired number of dimensions (default 4)

    Returns:
    X_adjusted: 4D xarray DataArray with dims (t, z, y, x)
    orig_dim: int, original number of dimensions of X

    '''
    cdim = 0
    orig_dim = X.ndim
    while X.ndim < n:
        new_dim = f'j{cdim}'
        X = X.expand_dims({new_dim: 1})
        cdim = cdim + 1
    return X, orig_dim
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

def expand_array(X: xr.DataArray, pad_op:str="edge", ndeg:int=1):
    '''
    Args:
    X: xarray DataArray with dims (t, z, y, x)
    pad_op: string, 'edge' or 'interp'.
    ndeg: int, number of cells to join in i, j

    Uses xarray.DataArray.pad
    pad array with linear extrapolation (pad_op='interp')
    or edge vaslues (pad_op='edge')
    to have the sizes y, x multiples of ndeg.

    also add a rim of size ndeg, so I'm sure to have
    land points at the N, S, E boundaries
    when degrading resolution
    Returns:
    X_padded: xarray DataArray with dims (t, z, y_padded, x_padded)
    '''
    jpt, jpk, jpj, jpi = X.shape
    # compute padding sizes, integers 0 <= d <= ndeg -1
    dj = (ndeg - jpj % ndeg) % ndeg
    di = (ndeg - jpi % ndeg) % ndeg

    padding_south = dj + ndeg # applying delta + the future border of one cell
    padding_north = ndeg      # just the border
    padding_east = 0          # nothing on Atlantic buffer
    padding_west = di + ndeg  # applying delta + the future border of one cell
    
    # generate pad_width dict
    pw = {'y':(padding_south, padding_north), 'x':(padding_east, padding_west)}

    if pad_op=='edge':
        X_padded = X.pad(pad_width = pw, mode='edge')
    elif pad_op=='interp':
        X_padded = X.pad(pad_width = pw, mode='constant', constant_values=np.nan)
        X_padded = fill_rim_lin_2D(X_padded)
    return X_padded

def xpnd_wrap(X: xr.DataArray, pad_op:str="edge", ndeg:int=1):
    '''
    Args:
    X: xarray DataArray with dims (... y, x) or less
    pad_op: string, 'edge' or 'interp'.
    ndeg: int, number of cells to join in i, j

    Returns:
    X_expanded: xarray DataArray with dims (... y_padded, x_padded)
    '''
    X_4d, orig_dims = adjust_dims(X)
    X_4d_expand = expand_array(X_4d, pad_op, ndeg)
    X_expand = deadjust_dims(X_4d_expand, orig_dims)
    return X_expand

def e3_2_gdep(M, tuvw):
    '''
    reconstructs gdept from e3t_0
    cause no gdept in the original Nemo mesh
    '''
    e3 = M[f'e3{tuvw}_0']
    gdep = e3.cumsum(dim="z") - (e3 / 2.0)
    return gdep


# DEGRADE RESOLUTION





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

def degr_wrap(X24:xr.DataArray, degr_op:callable, ndeg=1, W=1.0):
    '''
        Wrapper to degrade resolution of 2D/3D/4D arrays.
        
        Parameters
        ----------
        X24 : xarray.DataArray
            Input array with dimensions (..., y, x).
        degr_op : Callable
            Regridding function that performs the degradation operation.
            Expected signature: degr_op(X24_dims, ndeg, W) -> xarray.DataArray
        ndeg : int, optional
            Number of cells to join in i, j directions. Default is 1.
        W : xarray.DataArray, optional
            Weights for weighted means computation. Default is 1.0.
        
        Returns
        -------
        xarray.DataArray
            Degraded array with the same dimensions as the input X24.
        
        Notes
        -----
        The function internally transforms the input to 4D format, applies
        the degradation operation on reshaped blocks, and then restores
        the original dimensionality.    

    '''
    X24_4D, nd = adjust_dims(X24) # make 4D
    X24_6D = reshape_blocks(X24_4D, ndeg)
    X_degr = degr_op(X24_6D, ndeg, W)
    X_degr = deadjust_dims(X_degr, nd) # back to original dims
    return X_degr

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


def degrade_mesh(M1:xr.DataArray, thresh:int=1, ndeg:int=1):
    '''
    generates new reduced mesh for each variable

    Arguments:
    M1 : xr.DataArray where for example e2t field is
        M1["e2t"] = xpnd_wrap(M["e2t"], 'edge', ndeg)
        so, it has dimensions (... y,x ) multipes of ndeg
    thresh : integer, number of water points needed to classify a cell as water
    ndeg : integer, number of cells to join in i, j directions

    Returns:
    M2: xr.DataArray
        with dimensions y/ndeg, x/ndeg
    '''
    M2 = {}
    #
    M2["e1f"] = degr_wrap(M1["e1f"], rg.isum_jstep, ndeg, W=None)
    M2["e1t"] = degr_wrap(M1["e1t"], rg.isum_jmean, ndeg, W=None)
    M2["e1u"] = degr_wrap(M1["e1u"], rg.isum_jmean, ndeg, W=None)
    M2["e1v"] = degr_wrap(M1["e1v"], rg.isum_jstep, ndeg, W=None)
    M2["e2f"] = degr_wrap(M1["e2f"], rg.jsum_istep, ndeg, W=None)
    M2["e2t"] = degr_wrap(M1["e2t"], rg.jsum_istep, ndeg, W=None)
    M2["e2u"] = degr_wrap(M1["e2u"], rg.jsum_istep, ndeg, W=None)
    M2["e2v"] = degr_wrap(M1["e2v"], rg.jsum_istep, ndeg, W=None)
    M2["e3t_0"] = degr_wrap(M1["e3t_0"], vwmean, ndeg, W=M1['V'])
    M2["e3u_0"] = degr_wrap(M1["e3u_0"], uawmean_istep, ndeg, W=M1['Au'])
    M2["e3v_0"] = degr_wrap(M1["e3v_0"], vawmean_jstep, ndeg, W=M1['Av'])
    M2["e3w_0"] = degr_wrap(M1["e3w_0"], vwmean, ndeg, W=M1['V'])

    M2["gdept"] = degr_wrap(M1["gdept"], vwmean, ndeg, W=M1['V'])
    M2["gdepw"] = degr_wrap(M1["gdepw"], vwmean, ndeg, W=M1['V'])
    M2["glamf"] = degr_wrap(M1["glamf"], rg.jmean_istep, ndeg, W=None)
    M2["glamt"] = degr_wrap(M1["glamt"], rg.imean_jstep, ndeg, W=None)
    M2["glamu"] = degr_wrap(M1["glamu"], rg.jmean_istep, ndeg, W=None)
    M2["glamv"] = degr_wrap(M1["glamv"], rg.imean_jstep, ndeg, W=None)
    M2["gphif"] = degr_wrap(M1["gphif"], rg.imean_jstep, ndeg, W=None)
    M2["gphit"] = degr_wrap(M1["gphit"], rg.jmean_istep, ndeg, W=None)
    M2["gphiu"] = degr_wrap(M1["gphiu"], rg.jmean_istep, ndeg, W=None)
    M2["gphiv"] = degr_wrap(M1["gphiv"], rg.imean_jstep, ndeg, W=None)
    M2["nav_lat"] = degr_wrap(M1["nav_lat"], rg.jmean_istep, ndeg, W=None)
    M2["nav_lon"] = degr_wrap(M1["nav_lon"], rg.imean_jstep, ndeg, W=None)
    M2["nav_lev"] = M1["nav_lev"]
    M2["tmask"] = waterpt_thresh(M1["tmask"], thresh ,ndeg)
    M2['umask'] = degr_wrap(M1['umask'], rg.jany_istep, ndeg, W=None)
    M2['vmask'] = degr_wrap(M1['vmask'], rg.iany_jstep, ndeg, W=None)

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
    # zero at x=0 (western rim)
    MMed['tmask'].loc[dict(x=0)] = 0
    print('end ZW')
    return MMed


if __name__=='__main__':
    data = np.array([[1, 2],
                     [3, 4],
                     [5, 6]])   # shape (3, 2)
    arr= xr.DataArray(data, dims=('x','y'), coords={'x':[0,1,2], 'y':[10,20]})
    arr.pad({'y':(2,2),'x':(2,2)},mode='edge')
    X, nd = adjust_dims(arr)
    print(X)

    X_exp = xpnd_wrap(arr, pad_op='edge', ndeg=3) # su array piccoli fa una cosa strana
    data = np.ones((20,20), dtype=np.float32)
    arr= xr.DataArray(data, dims=('x','y'), coords={'x':np.arange(20), 'y':np.arange(20,40)})
    print(xpnd_wrap(arr, pad_op='edge', ndeg=4).sizes)

    X = xpnd_wrap(arr, pad_op='edge', ndeg=3)
    X_4D, nd = adjust_dims(X)
    X24_6D = reshape_blocks(X_4D, ndeg=3)
    X_degr = rg.jsum_istep(X24_6D, ndeg=3, W=None)
    if False:
        M1=load_mesh('meshmask.nc', ndeg=6)
        umask6d = reshape_blocks(M1["umask"], 6)



