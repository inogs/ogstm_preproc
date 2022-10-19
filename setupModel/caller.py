from commons.genUserDateList import getTimeList, getCouplesTimeList
import writeSerialJob

# USER SECTION #####################################
AVE1 = getTimeList("20190101-00:00:00" , "20230101-00:00:00", days=1 )    

AVE2 = getTimeList("20190107-00:00:00" , "20230101-00:00:00" ,days=7 )

RST = getTimeList("20190107-00:00:00"  , "20230101-12:00:00", days=7  )

###################################################


dateFormat="%Y%m%d-%H:%M:%S"

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


import sys
sys.exit()
# This section is not maintained

TL  = getCouplesTimeList("20130101-00:00:00"     , "20140101-00:00:00",months =2 , seconds=0 )
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
