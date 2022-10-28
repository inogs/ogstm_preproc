#! /bin/bash

#SBATCH --job-name=ATM
#SBATCH -N1
#SBATCH --ntasks-per-node=10
#SBATCH --time=01:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS21_PRACE_P
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg

cd $SLURM_SUBMIT_DIR
. ../../profile.inc
source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh 
unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false

INPUTFILE=/g100_work/OGS_prodC/OPA/V9C-dev/inpdir/forecast/20221027/ECMWF/20221027_20221106-ECMWF---AM0100-MEDATL-b20221027_fc00-fv11.00.nc 
OUTDIR=/g100_scratch/userexternal/gbolzon0/V9C/ONLINE_ATM/out 
MASKFILE=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc

mkdir -p $OUTDIR
my_prex_or_die "mpirun -np 10 python forecast_input_gen.py -i $INPUTFILE -o $OUTDIR -m $MASKFILE"
