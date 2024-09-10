import numpy as np
from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
import datetime
from bitsea.commons.timerequestors import Clim_day

import bitsea.commons.genUserDateList as DL
daily=DL.getTimeList("19970601-00:00:00", "20010502-12:00:00", "days=1")
TL = TimeList(daily)

deltadays=15
basetime=datetime.datetime(2000,1,1)
for julian in range(0,1):
    central = basetime + datetime.timedelta(days=julian)
    
    i_start = julian - deltadays
    i_end   = julian + deltadays
    index_list = []
    for iday in range(i_start,i_end):
        day = basetime + datetime.timedelta(days=iday)
        #print day
        clim_day_req = Clim_day(day.month,day.day)
        ii,w = TL.select(clim_day_req)
        index_list.extend(ii)  
    indsort=np.argsort(index_list)
    II = [index_list[k] for k in indsort]
    timelist= [TL.Timelist[k] for k in II]
    nFrames = len(II)
    for k in range(nFrames):
        print II[k],timelist[k]
        
    #starttime = basetime + datetime.timedelta(days = julian - deltadays)
    #end__time = basetime + datetime.timedelta(days = julian + deltadays)
    #print starttime,end__time