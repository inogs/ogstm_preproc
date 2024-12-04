#! /bin/bash

#SBATCH --job-name=Po
#SBATCH -N2
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_devC
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh
unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false
source ../profile.inc

MASKFILE=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc
INPUTBC=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/BC/inputs

python write_VPnutrients.py -i $INPUTBC -m $MASKFILE
python main.py

RESTARTSDIR=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/IC/RST_2018
OUTDIR=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/BC/out

cd atlantic
python atlantic_generator.py -i $INPUTBC -o $OUTDIR -m $MASKFILE --rst $RESTARTSDIR

cd ..
CMCC_MASK=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask_CMCC.nc
FORCINGS=/g100_scratch/userexternal/gbolzon0/V9C/2019/FORCINGS/
my_prex_or_die "mpirun python Po/online/po_generation.py -i $FORCINGS -o $OUTDIR -M $CMCC_MASK -m $MASKFILE"

