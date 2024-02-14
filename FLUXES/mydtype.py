#! /usr/bin/python

import numpy as np


trans_dt=np.dtype([('name','S12'),('type','S12'),('lon1',np.float32),('lon2',np.float32),\
                                                ('lat1',np.float32),('lat2',np.float32),\
                                                ('depth1',np.float32),('depth2',np.float32)])


ind2d   =np.dtype([('lat',int),('lon',int)])


# datatype for fluxes

flux_dt =np.dtype([('adv-u',np.float32),('adv-v',np.float32),('adv-w',np.float32),('sed-w',np.float32),\
                   ('hdf-x',np.float32),('hdf-y',np.float32),('zdf-z',np.float32)])

