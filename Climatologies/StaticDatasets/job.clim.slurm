#!/bin/bash

#SBATCH --job-name=static
#SBATCH -N1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=100gb
#SBATCH --account=OGS23_PRACE_IT 
#SBATCH --partition=g100_usr_prod

cd $SLURM_SUBMIT_DIR
. ../../profile.inc
module purge
module load profile/base
module load autoload
module load intelmpi/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary netcdf/4.7.4--oneapi--2021.2.0-ifort ncview
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18_new/bin/activate
#export PYTHONPATH=/g100_scratch/userexternal/vdibiagi/bit.sea

my_prex_or_die "python clim_from_static_open.py"

MASKFILE=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc
OUTDIR=/g100_work/OGS_test2528/Observations/TIME_RAW_DATA/STATIC/MedBGCins/Clim_Annual_Ref/
my_prex_or_die "python clim_for_medbgcins_archive.py -m $MASKFILE -o $OUTDIR"

