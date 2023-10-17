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
PERCENTILES_DIR=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/percentiles/
         PUBDIR=/g100_work/OGS_devC/Benchmark/pub/Benchmark/eas7_V10C/PHYS

RUNNAME=benchmark

mkdir -p $METRICS_2D $PERCENTILES_DIR $PUBDIR
my_prex_or_die "mpirun python metrics_2d.py -i $INPUTDIR -o $METRICS_2D -m $INGVMASK"

my_prex_or_die "python Kz_anom_synthesis.py -i $METRICS_2D -m $INGVMASK -o $PUBDIR"
my_prex_or_die "python metrics_2d_percentiles.py -i $METRICS_2D -m $INGVMASK -o $PERCENTILES_DIR"
my_prex_or_die "python metrics_2d_percentiles_plot.py -i $PERCENTILES_DIR -o $PUBDIR -n $RUNNAME"

# my_prex "python metrics_2d_percentiles_plot_multirun.py" # to plot in https://medeaf.ogs.it/internal-validation/Benchmark/eas7_V10C/multirun_vs_eas8/

OUTPUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/FORCINGS/DELTAT
my_prex_or_die "mpirun python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK"

