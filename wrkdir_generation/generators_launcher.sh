#! /bin/bash

HERE=$PWD
MONTH=11

month=`printf "%02d" $MONTH`

for TYPE in T0 T1A  T1B1 ; do
   DIR=R${month}/${TYPE}/wrkdir/MODEL
   echo $DIR 
   mkdir -p $DIR
   cp generator.sh $DIR
   cd $DIR
   ./generator.sh -m $MONTH -f /gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/FORCINGS/READY_FOR_MODEL/${TYPE}
   
   sbatch --job-name=R02/${TYPE} job.slurm 
   cd $HERE
done



