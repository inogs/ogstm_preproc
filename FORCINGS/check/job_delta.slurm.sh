#! /bin/bash

#SBATCH --job-name=deltat
#SBATCH -N2
#SBATCH --ntasks-per-node=30
#SBATCH --time=02:00:00
#SBATCH --mem=300gb
#SBATCH --account=OGS21_PRACE_P
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg

cd $SLURM_SUBMIT_DIR
. ../../profile.inc

source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh

unset I_MPI_PMI_LIBRARY


INPUTDIR=/g100_scratch/userexternal/gbolzon0/V9C/2019/FORCINGS/
OUTPUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/FORCINGS/DELTAT
INGVMASK=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask_CMCC.nc


my_prex_or_die "mpirun python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK"

