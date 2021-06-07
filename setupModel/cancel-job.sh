#!/bin/bash

##################################################
#           cancel-job.sh 
#    Cancels all the job with jobid > argument_1
# 
#    Author : Paolo Lazzari
# ################################################

if [ -z "$1" ] ; then
    echo "Minimum Job Number argument is required.  Run as '$0 jobnum'"
    exit 1
fi

minjobnum="$1"

for j in $(squeue --user="$USER" --noheader --format=%i) ; do
  if [ "$j" -gt "$minjobnum" ] ; then
    scancel "$j"
  fi
done
