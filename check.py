import config as conf
from bclib.river import river
from bclib.bounmask import bounmask

BOUN = bounmask(conf)
index_inv = BOUN.load('index_inv')



filename=conf.dir_out + "TIN_20040515-00:00:00.nc"
R=river(conf)
idxt = R.load_from_file(filename, 'riv_idxt')


for p in idxt:
    print index_inv[p,:]
