import os

JOBS=[]
filein=file("SubmissionScheme","r")
for line in filein:
    JOBS.append(line[:-1])
filein.close()   


def writeSerialJob(NSTEPS_LAUNCHING_MODEL):
    
    deflines=[]
    
    deflines.append("#! /bin/bash " + "\n"*2)
    deflines.append("DEPENDENCE=''" + "\n"*2)
    
    for k in xrange(NSTEPS_LAUNCHING_MODEL):
        for jobdefline in JOBS:            
            line="JOBID=`" + jobdefline + "` "
            deflines.append(line%(k,k,k));
            line="DEPENDENCE=\" --dependency=afterany:$JOBID \" " 
            deflines.append(line)
            
        #line= "JOBID=`qsub $DEPENDENCE -o job.out.%d -e job.err.%d jobS.pbs` " %(k,k)
        
        
        

        
        
        
    
    
    outfilename="launcher.sh"    
    outfile=file(outfilename,"w")
    for line in deflines:
        outfile.writelines(line + "\n")
    outfile.close()
    os.system("chmod 744 "+outfilename)


