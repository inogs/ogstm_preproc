#!/bin/bash

#SBATCH --job-name=BC
#SBATCH -N2
#SBATCH --ntasks-per-node=32
#SBATCH --time=0:10:00
#SBATCH --mem=300gb
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_dbg


cd $SLURM_SUBMIT_DIR
. ../../profile.inc

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh
unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false

INPUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FORCINGS/READY_FOR_MODEL/
OUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/BC/out

mkdir -p $OUTDIR

INGVMASK=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc
MASKFILE=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc

my_prex_or_die "mpirun python river_generation.py -i $INPUTDIR -o $OUTDIR -m $MASKFILE -M $INGVMASK"

