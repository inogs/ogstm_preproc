from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import Sat.SatManager as Sat
import numpy as np

INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"
TI = TimeInterval("19500101","20500101","%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')


COUNT = np.zeros((1580,3308),np.int32)
for ifile, filename in enumerate(TL.filelist):
    print ifile
    fillvalue=-999
    A=Sat.readfromfile(filename,'KD490')
    ii=np.isnan(A)
    A[ii] = fillvalue
    filled = A==fillvalue
    goods = ~filled 
    COUNT[goods] = COUNT[goods]+1
    


np.save('Kd_map_occurency',COUNT)