#! /bin/bash
  
#SBATCH --job-name=hr_avg
#SBATCH -N1
#SBATCH --ntasks-per-node=48
#SBATCH --time=1:00:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_test2528
#SBATCH --partition=g100_usr_prod

cd $SLURM_SUBMIT_DIR
. ../profile.inc

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh

unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false

MASKFILE=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc
INPUTDIR=/g100_work/OGS_test2528/NECCTON/OPTICS/2022
OUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/OPTICS/READY_FOR_MODEL/
my_prex_or_die "mpirun -np 48 python hourly_averager.py -m $MASKFILE -i $INPUTDIR -o $OUTDIR"

INPUTDIR=/g100_work/OGS_test2528/NECCTON/OPTICS/2023
OUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/OPTICS/READY_FOR_MODEL/
my_prex_or_die "mpirun -np 48 python hourly_averager.py -m $MASKFILE -i $INPUTDIR -o $OUTDIR"
