#! /bin/bash

####    tagliafile.sh              #           #####
#   cut a single monthly file in daily files       #
#                                                  #
#   Author: GB. 2019.02.11                         #
####################################################


usage() {
echo "Generates daily files from a single monthly file "
echo "SYNOPSYS"
echo "tagliafile.sh [ -i INPUTDIR] [-o outputdir ]"
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
      "-o" ) DAILY_DIR=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done


filename=`basename $inputfile`
var=${filename:36:1}
ndays=${filename:28:2}
YEAR=${filename:22:4}
MONTH=${filename:26:2}


for (( day=1; day < ndays+1; day++ )) ; do
          DAY=`printf %02d $day`
          DATE=${YEAR}${MONTH}${DAY}
          DAILY=${DAILY_DIR}/${DATE}

          ncks -O -F -d time_counter,$DAY,$DAY $inputfile  ${DAILY}${var}.nc
done
