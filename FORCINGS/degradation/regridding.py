import xarray as xr

# REGRIDDING FUNCTIONS
def isum_jmean(X: xr.DataArray, ndeg=1, W=1.0):
    '''
    sum along longitude (i) & mean along latitude (j)
    no need to weight by area! for e1t

    The algorithm on a single block works as follows:
    [ a11 a12 a13     --> sum these ones
      a21 a22 a23
      a31 a32 a33 ]

    For e1t, we do prefer to repeat the sum for each j_b, and then average along j_b

    Arguments:
    X: xarray DataArray with dims 
          (..., y,      y_b,  x,       x_b) with dims sizes
    
    Returns:
    Xr: xarray DataArray with dims (..., y,  x )

    '''
    X_yb = X.sum(dim='x_b', skipna=True)  #blocks sum along lon (i), returns (..., y, y_b, x)
    Xr = X_yb.mean(dim='y_b', skipna=True) #block mean along lat (j)
    return Xr

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

if __name__ == "__main__":
    import numpy as np
    
    data = np.ones((20,4,20,4), dtype=np.float32)
    arr= xr.DataArray(data, dims=('y','y_b','x','x_b'), coords={'x':np.arange(20), 'y':np.arange(20)})
    a=isum_jmean(arr, ndeg=4, W=None)

    data= np.arange(20).reshape((5,4)).astype(np.float32)
    arr= xr.DataArray(data, dims=('y_b','x_b'))
    isum_jmean(arr, ndeg=4, W=None)

                    
