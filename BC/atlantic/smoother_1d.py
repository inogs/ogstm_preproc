import numpy as np

def smoother_1d(y, x):
    ''' Simple smoother
    '''
    xNew=0.5*(x[:-1]+x[1:])
    yNew=np.interp(xNew,x,y)
    y=np.interp(x,xNew,yNew)
    return y


if __name__ == '__main__':
    x=np.arange(10)
    y=[0,0,0,0,0,4,4,4,4,4]
    print("y before interpolation: ",y)
    for itime in range(2):
        y=smoother_1d(y,x)
    print("and after: ", y)
