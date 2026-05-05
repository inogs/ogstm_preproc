#! /bin/bash

#SBATCH --job-name=ppcon

#SBATCH -N1
#SBATCH --ntasks-per-node=10
#SBATCH --time=1:00:00
#SBATCH --mem=50gb
#SBATCH --account=OGS_test2528
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo


. ../profile.inc

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh

BITSEA=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FLOATS/bit.sea
export PYTHONPATH=$BITSEA/src
cd $BITSEA/src/bitsea/Float/ppcon/

TRAIN_DIR=/g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/FLOAT/PPCON_TRAIN

export ONLINE_REPO=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FLOATS/ONLINE

ONLINE_REPO_CLUSTERING=${ONLINE_REPO}/PPCON/clustering/

mkdir -p $ONLINE_REPO_CLUSTERING


for YEAR in $( seq 2012 2026 ) ; do

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

cd ..

my_prex_or_die "mpirun python dump_index.py -i ${ONLINE_REPO}/PPCON -o ${ONLINE_REPO}/PPCON/Float_Index.txt -t ppcon_float"

