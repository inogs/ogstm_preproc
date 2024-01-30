#! /bin/bash

####    cut_singlefile_op.sh               #########
#   cut a single daily file in 6-hourly files      #
#    following eas8 operational names              #
#                                                  #
#   Author: GB. 2021.06.10                         #
####################################################


usage() {
echo "Generates 6-hourly files from a single daily file "
echo "SYNOPSYS"
echo "cut_singlefile_op.sh [ -i inputfile] [-o outputdir ]"
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
var_date8=${filename:0:9}

for (( iframe=1; iframe < 5; iframe++ )) ; do
	  hr=$((  $(( $iframe * 6 )) - 3 ))
	  HR=`printf %02d $hr`
          outputfile=${HOURLY_DIR}/${var_date8}-${HR}:00:00.nc

          ncks -O -F -d time_counter,$iframe,$iframe $inputfile $outputfile
done

