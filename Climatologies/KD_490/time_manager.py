import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import datetime
from commons.timerequestors import Clim_day
INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"
TI = TimeInterval("19500101","20500101","%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')

deltadays=15
basetime=datetime.datetime(2000,1,1)

def getfilelist(julian):
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

