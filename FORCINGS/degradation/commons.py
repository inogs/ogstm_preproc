import xarray as xr


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

def reshape_blocks(X, ndeg=1):
    '''
    reshape in blocks of shape ndeg along j, i
    xarray's convoluted way of doing np.reshape, esticatsi!
    assume I'm always reshaping dims called y, x
    chunking is necessary to avoid OOM
    (still cuz xarray's way of reshaping is not super smart)

    Args:
    X: xarray DataArray with dims (... y, x)
    ndeg: int, number of cells to join in y, x directions
          x and y sizes must be multiples of ndeg
    Returns:
    X_reshaped: xarray DataArray with dims 
          (..., y,      y_b,  x,       x_b) with dims sizes
          (..., y/ndeg, ndeg, x/ndeg, ndeg)
    So, X_reshaped[... 0,:,0,:] is the first block of size ndeg x ndeg

    '''
    cs = ndeg * 5
    X = X.chunk({'y':cs, 'x':cs})
    X_reshaped = X.coarsen(y=ndeg, x=ndeg).construct(y=('y', 'y_b'), x=('x', 'x_b'))

    return X_reshaped

def degrade_wrap(X, degr_op, B=1.0, ndeg=1):
    print('---1---')
    X = reshape_blocks(X, ndeg)
    print('---2---')
    X = degr_op(X, ndeg, B)
    print('...')
    return X

def dump_netcdf(M2: xr.Dataset, outfile: str):
    '''
    Generic function to save a xarray Dataset to a NetCDF file.
    Args:
    M2: xarray Dataset to be saved.
    outfile: str, path to the output NetCDF file.
    '''
    print(f'saving: {outfile}')
    M2.to_netcdf(outfile, mode='w', unlimited_dims='time', format='NETCDF4_CLASSIC')
    return