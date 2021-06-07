# The first part produces 1.aveTimes, 2.aveTimes and restartTimes, sequences of 17-chars dates

# CONDITION FOR DATA ASSIMILATION: 1.aveTimes list must contain daTimes, with hour set to 00:00
# E.g. if daTime = 20070712-12:00:00 ==>  1.aveTimes must contain 20070712-00:00:00
AVE_START_TIME_1="20130101-00:00:00"
AVE_END___TIME_1="20140101-00:00:00"
AVE_DELTA_1     ="days  =  1 "


AVE_START_TIME_2="20130101-00:00:00"
AVE_END___TIME_2="20140101-00:00:00"
AVE_DELTA_2     ="months=  1 "


RST_START_TIME="20130101-00:00:00"
RST_END___TIME="20140101-12:00:00"
RST_DELTA     ="months=1 "#"months = 1 "


# This 2nd part produces the Time Table JobTimeDivision.dat

TTb_START_TIME="20130101-00:00:00"
TTb_END___TIME="20140101-00:00:00"
TTb_DELTA     ="months  = 2 "
TTb_OVERLAP   ="seconds = 0"

# to avoid errors, the first column of time table (Start Times)
# must match with restartTimes

