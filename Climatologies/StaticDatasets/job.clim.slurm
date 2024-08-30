#!/bin/bash

#SBATCH --job-name=N1p
#SBATCH -N1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:49:00
#SBATCH --mem=100gb
#SBATCH --account=OGS23_PRACE_IT 
#SBATCH --partition=g100_usr_prod
#SBATCH --error N1p.err
#SBATCH --output N1p.out

cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load autoload
module load intelmpi/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary netcdf/4.7.4--oneapi--2021.2.0-ifort ncview
LIBDIR=/g100_work/OGS23_PRACE_IT/COPERNICUS/V10C/HOST/g100
export    HDF5_DIR=$LIBDIR
export NETCDF4_DIR=$LIBDIR
export    GEOS_DIR=$LIBDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBDIR

source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18/bin/activate
export PYTHONPATH=$PYTHONPATH:/g100_scratch/userexternal/vdibiagi/bit.sea

date

python /g100_scratch/userexternal/vdibiagi/preproc/Climatologies/StaticDatasets/clim_from_static_open.py 

date
exit 0



