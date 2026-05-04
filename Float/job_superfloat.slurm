#! /bin/bash

#SBATCH --job-name=Floats
#SBATCH -N1
#SBATCH --ntasks-per-node=10
#SBATCH --time=2:00:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_test2528
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

cd $SLURM_SUBMIT_DIR
. ../profile.inc


source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh

BITSEA=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FLOATS/bit.sea
export PYTHONPATH=$BITSEA/src

cd $BITSEA/src/bitsea/Float

export ONLINE_REPO=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FLOATS/ONLINE/
SUPERFLOAT_DIR=$ONLINE_REPO/SUPERFLOAT/
mkdir -p $SUPERFLOAT_DIR
COMMONS="-s 20120101 -e 20270101 -o $SUPERFLOAT_DIR"

my_prex_or_die "python superfloat_chla.py $COMMONS"
my_prex_or_die "cp -r $SUPERFLOAT_DIR /g100_scratch/userexternal/gbolzon0/V13C/SUPERFLOAT_CHLA/"

export CLIM_MED_BGC_INS_FILE=/g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/FLOAT/BGCINS/Climatology_N3n_600_800m.nc
my_prex_or_die "python superfloat_nitrate.py $COMMONS"
my_prex_or_die "cp -r $SUPERFLOAT_DIR /g100_scratch/userexternal/gbolzon0/V13C/SUPERFLOAT_NITRATE/"

mkdir -p $ONLINE_REPO/SUPERFLOAT/Oxydiag
my_prex_or_die "mpirun -np 10 python superfloat_oxygen.py $COMMONS -O $ONLINE_REPO/SUPERFLOAT/Oxydiag/"

my_prex_or_die "cp -r $SUPERFLOAT_DIR /g100_scratch/userexternal/gbolzon0/V13C/SUPERFLOAT_DOXY/"
my_prex_or_die "python superfloat_bbp700.py $COMMONS -f "
my_prex_or_die "python superfloat_kd490.py $COMMONS"
my_prex_or_die "python superfloat_par.py $COMMONS"
my_prex_or_die "python superfloat_ph.py $COMMONS"

my_prex_or_die "mpirun -np 10 python dump_index.py -i $SUPERFLOAT_DIR -o $SUPERFLOAT_DIR/Float_Index.txt -t superfloat"

