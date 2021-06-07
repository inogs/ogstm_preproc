# written to updated Start_End_Times for multiples runs of opa.xx
# submitted by launcher.sh 
import glob

TimeT = []
filein=file("JobTimeDivision.dat","r")
for line in filein:
     coppia = [line[0:17], line[18:-1] ]
     TimeT.append(coppia)
filein.close()


for inizio,fine in TimeT:
    print "searching for RST." + inizio + ".N1p.nc ........"   
    a=glob.glob("RESTARTS/RST." + inizio + ".N1p.nc")
    
    if a.__len__()==0:#exit condition
        break    
    Inizio=inizio
    Fine = fine
# then they are the final values
try:
    print Inizio
    print Fine
    
    
    StartEnd = file("Start_End_Times","w")
    StartEnd.write(Inizio + "\n")
    StartEnd.write(Fine   + "\n")
    StartEnd.close()
    print "Start_End_Times written"
except:

    print "No files have been found to match the requests of JobTimeDivision.dat"
    print " ******* Start_End_Times not written ********** "
    print "Please re-generate date lists or check RST files"
    print "writeStartEnd.py search for the beginning of the list, then ahead"





