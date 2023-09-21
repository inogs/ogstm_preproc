#! /bin/bash

#SBATCH --job-name=deltat
#SBATCH -N2
#SBATCH --ntasks-per-node=30
#SBATCH --time=02:00:00
#SBATCH --mem=300gb
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg

cd $SLURM_SUBMIT_DIR
. ../../profile.inc

source /g100_work/OGS21_PRACE_P/COPERNICUS/sequence3.sh

unset I_MPI_PMI_LIBRARY


INPUTDIR=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/READY_FOR_MODEL
INGVMASK=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask_CMCC.nc


METRICS_2D=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/output
OUTDIR=/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/synthesis/
PERCENTILES_DIR=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/Ved_percentiles/
RUNNAME=benchmark
mkdir -p $METRICS_2D $OUTDIR $PERCENTILES_DIR
my_prex_or_die "mpirun python metrics_2d.py -i $INPUTDIR -o $METRICS_2D -m $INGVMASK"

my_prex_or_die "python anom_synthesis.py -i $METRICS_2D -m $INGVMASK -o $OUTDIR"
my_prex_or_die "python metrics_2d_percentiles.py -i $METRICS_2D -m $INGVMASK -o $PERCENTILES_DIR"
my_prex_or_die "python Kz_statplot.py -i $PERCENTILES_DIR -o $OUTDIR -n $BENCHMARK"
 

exit 0
OUTPUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/FORCINGS/DELTAT
my_prex_or_die "mpirun python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK"

