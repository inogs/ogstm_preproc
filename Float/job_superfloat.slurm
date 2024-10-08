#! /bin/bash

#SBATCH --job-name=Floats
#SBATCH -N1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=50gb
#SBATCH --account=OGS_devC
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

cd $SLURM_SUBMIT_DIR
. ../profile.inc


module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
module load netcdf/4.7.4--oneapi--2021.2.0-ifort
module load netcdff/4.5.3--oneapi--2021.2.0-ifort
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18_new/bin/activate

BITSEA=/g100_scratch/userexternal/gbolzon0/V11C/bit.sea
export PYTHONPATH=$BITSEA/src

cd $BITSEA/src/bitsea/Float

export ONLINE_REPO=/g100_work/OGS_devC/V11C/TRANSITION/ONLINE/
SUPERFLOAT_DIR=/g100_scratch/userexternal/gbolzon0/V11C/SUPERFLOAT/
COMMONS="-s 20120101 -e 20250101 -o $SUPERFLOAT_DIR"

my_prex_or_die "python superfloat_chla.py $COMMONS"
my_prex_or_die "cp -r $SUPERFLOAT_DIR /g100_scratch/userexternal/gbolzon0/V11C/SUPERFLOAT_CHLA/"

my_prex_or_die "python superfloat_nitrate.py $COMMONS"
my_prex_or_die "cp -r $SUPERFLOAT_DIR /g100_scratch/userexternal/gbolzon0/V11C/SUPERFLOAT_NITRATE/"

