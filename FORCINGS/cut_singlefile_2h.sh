#! /bin/bash

####    cut_singlefile_2h.sh               #########
#   cut a single daily file in bi-hourly files     #
#                                                  #
#   Author: GB. 2021.06.10                         #
####################################################


usage() {
echo "Generates bi-hourly files from a single daily file "
echo "SYNOPSYS"
echo "cut_singlefile_2h.sh [ -i inputfile] [-o outputdir ]"
echo "uses ncks"
echo ""
}


if [ $# -lt 4 ] ; then
  usage
  exit 1
fi

for I in 1 2 ; do
   case $1 in
      "-i" ) inputfile=$2;;
      "-o" ) HOURLY_DIR=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done


filename=`basename $inputfile`
var=${filename:${#filename}-4:1}
date8=${filename:${#filename}-18:8}

for (( iframe=1; iframe < 13; iframe++ )) ; do
	  hr=$((  $(( $iframe * 2 )) - 1 ))
	  HR=`printf %02d $hr`
          outputfile=${HOURLY_DIR}/${var}${date8}-${HR}:00:00.nc

          ncks -O -F -d time_counter,$iframe,$iframe $inputfile $outputfile
done
