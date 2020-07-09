import numpy as np
import datetime
from commons.timerequestors import Clim_day

deltadays=15
basetime=datetime.datetime(2001,1,1)

def getfilelist(julian, TL):
    #central = basetime + datetime.timedelta(days=julian)
    
    i_start = julian - deltadays
    i_end   = julian + deltadays
    index_list = []
    for iday in range(i_start,i_end):
        day = basetime + datetime.timedelta(days=iday)
        #print day
        clim_day_req = Clim_day(day.month,day.day)
        ii,_ = TL.select(clim_day_req)
        index_list.extend(ii)  
    indsort=np.argsort(index_list)
    II = [index_list[k] for k in indsort]
    filelist = [TL.filelist[k] for k in II]
    return II, filelist

