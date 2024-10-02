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

INPUTDIR=/g100_scratch/userexternal/gbolzon0/EFAS/FORCINGS/SIMU_EFAS_V5/READY_FOR_MODEL/
OUTDIR=/g100_work/OGS_devC/V11C/RIVERS/RUNS_SETUP/run06/out

mkdir -p $OUTDIR

INGVMASK=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask_CMCC.nc
MASKFILE=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc

my_prex_or_die "mpirun python river_generation.py -i $INPUTDIR -o $OUTDIR -m $MASKFILE -M $INGVMASK"

