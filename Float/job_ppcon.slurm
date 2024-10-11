#! /bin/bash

#SBATCH --job-name=ppcon
#SBATCH --ntasks=1
#SBATCH --time=02:05:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_all_serial


. ../profile.inc

module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
module load netcdf/4.7.4--oneapi--2021.2.0-ifort
module load netcdff/4.5.3--oneapi--2021.2.0-ifort
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18_new/bin/activate

BITSEA=/g100_scratch/userexternal/gbolzon0/V11C/bit.sea
cd $BITSEA/src/bitsea/Float/ppcon/

TRAIN_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/PPCON_TRAIN

export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V11C/
export ONLINE_REPO=/g100_scratch/userexternal/gbolzon0/V11C
ONLINE_REPO_CLUSTERING=${ONLINE_REPO}/PPCON/clustering/

mkdir -p $ONLINE_REPO_CLUSTERING


for YEAR in $( seq 2012 2024 ) ; do

    FILE_INPUT=${ONLINE_REPO}/PPCON/Float_Index_${YEAR}.txt
    grep ${YEAR} ${ONLINE_REPO}/SUPERFLOAT/Float_Index.txt > $FILE_INPUT

    # copying from superfloat
    for wmodir in $( cat ${FILE_INPUT} | cut -d "/" -f1 | uniq ) ; do
      mkdir -p ${ONLINE_REPO}/PPCON/${wmodir}
    done
    for filename in $( cat ${FILE_INPUT} | cut -d "," -f1 ) ; do    
       cp ${ONLINE_REPO}/SUPERFLOAT/$filename ${ONLINE_REPO}/PPCON/${filename}
    done



    my_prex_or_die "python -u clustering/clustering.py -i $ONLINE_REPO -u $FILE_INPUT -o $ONLINE_REPO_CLUSTERING"

    my_prex_or_die "python -u make_generated_ds/generate_netcdf_netcdf4.py -t $ONLINE_REPO_CLUSTERING -m $TRAIN_DIR -p $ONLINE_REPO/PPCON"


done


my_prex_or_die "python dump_index.py -i $ONLINE_REPO${PPCON} -o $ONLINE_REPO${PPCON}/Float_Index.txt -t ppcon_float"

