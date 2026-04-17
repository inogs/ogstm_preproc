#! /bin/bash

#SBATCH --job-name=floatClim
#SBATCH -N1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:20:00
#SBATCH --mem=300gb
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod


module load autoload
module load intelmpi/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary netcdf/4.7.4--oneapi--2021.2.0-ifort ncview
LIBDIR=/g100_work/OGS23_PRACE_IT/COPERNICUS/V10C/HOST/g100
export    HDF5_DIR=$LIBDIR
export NETCDF4_DIR=$LIBDIR
export    GEOS_DIR=$LIBDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBDIR

cd $SLURM_SUBMIT_DIR

. ../profile.inc

source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18_new/bin/activate
export ONLINE_REPO=/g100/home/userexternal/camadio0/CLIMATOLOGIE/FLOAT_CLIMATOLOGIES/new_clim_FLOAT_V12C_2012_2024/


export PYTHONPATH=$PYTHONPATH:/g100/home/userexternal/camadio0/CLIMATOLOGIE/FLOAT_CLIMATOLOGIES/new_clim_FLOAT_V12C_2012_2024/bit.sea/src/

VARNAME='O2o'
OUTDIR=/g100/home/userexternal/camadio0/CLIMATOLOGIE/FLOAT_CLIMATOLOGIES/new_clim_FLOAT_V12C_2012_2024/

mkdir -p $OUTDIR
# CA per runquid2025 executed in a subset of SUPERFLOAT DATASET 2012-2024 
my_prex "python Yr_Climfloat_netcdf.py -o $OUTDIR -v $VARNAME"
