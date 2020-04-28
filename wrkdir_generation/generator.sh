#! /bin/bash

usage() {
echo "Generates wrkdir/MODEL by linking "
echo "SYNOPSYS"
echo "generator.sh [ -f FORCINGS_DIR] [-m month ]"
echo "EXAMPLE"
echo 'generator.sh -f /gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/FORCINGS/READY_FOR_MODEL/T0B/ -m 2 '
echo ""
}

if [ $# -lt 4 ] ; then
  usage
  exit 1
fi

for I in 1 2 ; do
   case $1 in
      "-f" ) FORCINGS_DIR=$2;;
      "-m" ) MONTH=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done

ln -fs $FORCINGS_DIR FORCINGS

STARTTIME=`printf "2018%02d01-00:00:00" $MONTH `
END__TIME=`printf "2018%02d16-00:00:00" $MONTH `

echo -e "${STARTTIME}\n${END__TIME}" > Start_End_Times

mkdir -p AVE_FREQ_1 AVE_FREQ_2 RESTARTS

#-----------------------------------------------------------
cd RESTARTS

INITDIR=/gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/INIT
month=`printf "%02d" $MONTH ` 
for I in `ls $INITDIR/RST.2018${month}*nc `; do
  ln -fs $I
done

cd ..
#-----------------------------------------------------------
ln -fs /gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/CODE/OGSTM_BUILD/ogstm.xx

for I in `ls /gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/CODE/ogstm/ready_for_model_namelists/*nml `; do
   ln -fs $I
done

for I in `ls -d /gpfs/scratch/userexternal/gbolzon0/TESTS_0012/PREPROC/wrkdir/MODEL/*`; do
   ln -fs $I
done
# BC  bounmask.nc  domdec.txt  genInputsDatelists.sh  meshmask.nc  namelist.init  namelist.passivetrc

./genInputsDatelists.sh
