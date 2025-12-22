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
