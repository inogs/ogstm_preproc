#! /bin/bash

#SBATCH --job-name=co2
#SBATCH -N1
##SBATCH -N2
#SBATCH --ntasks-per-node=1
##SBATCH --ntasks-per-node=48
##SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --mem=30gb
#SBATCH --account=OGS_test2528_0
#SBATCH --partition=lrd_all_serial
##SBATCH --partition=dcgp_usr_prod
##SBATCH --partition=g100_usr_prod

# usage:
# copy this script to your directory,
# adjust to your own paths and HPG account
# submit with, sbatch main_co2.slurm

cd $SLURM_SUBMIT_DIR

##source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh
#source /leonardo_work/OGS23_PRACE_IT_0/COPERNICUS/sequence.sh

module load intel-oneapi-compilers/2023.2.1 netcdf-fortran/4.6.1--intel-oneapi-mpi--2021.10.0--oneapi--2023.2.0
# Initialize Conda (regrid_restart.py needs xesmf package)
source /leonardo/home/userexternal/ggalli00/miniconda3/etc/profile.d/conda.sh
conda activate /leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/PYENVS/regrid_env
export PYTHONPATH=/leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/git/bit.sea/src:$PYTHONPATH

gitdir=/leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/git/ogstm_preproc/
cd $gitdir/BC/bclib

unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false
source ../../profile.inc

my_prex_or_die "python co2.py"

