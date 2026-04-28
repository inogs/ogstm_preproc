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

MASKFILE=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc
INPUTBC=/g100_work/OGS_test2528/Benchmark/SETUP/PREPROC/BC/inputs/

python write_VPnutrients.py -i $INPUTBC -m $MASKFILE
python main.py

RESTARTSDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/IC/RST/
OUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/BC/out

cd atlantic
my_prex_or_die "python atlantic_generator.py -i $INPUTBC -o $OUTDIR -m $MASKFILE --rst $RESTARTSDIR"


cd ../EFAS
CMCC_MASK=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc
FORCINGS=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FORCINGS/READY_FOR_MODEL/
my_prex_or_die "mpirun python river_generation.py -i $INPUTDIR -o $OUTDIR -m $MASKFILE -M $CMCC_MASK"

