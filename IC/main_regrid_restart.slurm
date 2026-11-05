#!/bin/bash

#SBATCH --job-name=rst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem=24gb
#SBATCH --account=OGS23_PRACE_IT_0
##SBATCH --partition=dcgp_usr_prod
#SBATCH --partition=lrd_all_serial 

cd $SLURM_SUBMIT_DIR

module load intel-oneapi-compilers/2023.2.1 netcdf-fortran/4.6.1--intel-oneapi-mpi--2021.10.0--oneapi--2023.2.0

# Initialize Conda (regrid_restart.py needs xesmf package)
source /leonardo/home/userexternal/ggalli00/miniconda3/etc/profile.d/conda.sh
conda activate /leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/PYENVS/regrid_env

gitdir=/leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/git/ogstm_preproc/
cd $gitdir/IC

python regrid_restart.py -y regrid_restart.yaml
