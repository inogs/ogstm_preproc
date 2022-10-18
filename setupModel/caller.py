from commons.genUserDateList import *
import writeSerialJob
from user import *

dateFormat="%Y%m%d-%H:%M:%S"
AVE1 = getTimeList(       AVE_START_TIME_1      , AVE_END___TIME_1,days=1  )    
AVE2 = getTimeList(       AVE_START_TIME_2      , AVE_END___TIME_2,days=7 )

RST = getTimeList(       RST_START_TIME      , RST_END___TIME,days=7  )
TL  = getCouplesTimeList(TTb_START_TIME      , TTb_END___TIME,months =2 , seconds=0 )



outfile=open("1.aveTimes","w")
for dumptime in AVE1:
    outfile.write(dumptime.strftime(dateFormat) + "\n")
outfile.close()
print("1.aveTimes written")

outfile=open("2.aveTimes","w")
for dumptime in AVE2:
    outfile.write(dumptime.strftime(dateFormat) + "\n")
outfile.close()
print("2.aveTimes written")



outfile=open("restartTimes","w")
for dumptime in RST:
    outfile.write(dumptime.strftime(dateFormat) + "\n")
outfile.close()
print("restartTimes written")


NSTEPS_LAUNCHING_MODEL = len(TL)

print('NSTEPS_LAUNCHING_MODEL=', NSTEPS_LAUNCHING_MODEL)
writeSerialJob.writeSerialJob(NSTEPS_LAUNCHING_MODEL)
print("launcher.sh written")


outfile=open("JobTimeDivision.dat","w")
for tstart, tend__ in TL: 
    outfile.write(tstart.strftime(dateFormat) + " " + tend__.strftime(dateFormat) + "\n")
    #print tstart.strftime(dateFormat), "-->", tend__.strftime(dateFormat)
outfile.close()
print("JobTimeDivision.dat written")
