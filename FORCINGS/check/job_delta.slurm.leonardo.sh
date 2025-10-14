#! /bin/bash

#SBATCH --job-name=deltat
#SBATCH -N1
#SBATCH --ntasks-per-node=60
#SBATCH --time=00:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_test2528_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=dcgp_qos_dbg

cd $SLURM_SUBMIT_DIR
. ../../profile.inc

source /leonardo_work/OGS23_PRACE_IT_0/COPERNICUS/sequence.sh

unset I_MPI_PMI_LIBRARY


INPUTDIR=/leonardo_work/OGS_test2528_0/COPERNICUS/V12C/Forcings_Mixing/SIMU8/READY_FOR_MODEL/
INGVMASK=/leonardo_work/OGS_prod2528_0/OPA/V11C-prod/etc/static-data/MED24_141/meshmask.nc


     METRICS_2D=/leonardo_work/OGS_test2528_0/COPERNICUS/V12C/Forcings_Mixing/metrics/SIMU8/output
PERCENTILES_DIR=/leonardo_work/OGS_test2528_0/COPERNICUS/V12C/Forcings_Mixing/metrics/SIMU8/percentiles
         PUBDIR=/leonardo_work/OGS_test2528_0/pub/lfeudale/Forcings_Mixing/SIMU8

RUNNAME=benchmark

mkdir -p $METRICS_2D $PERCENTILES_DIR $PUBDIR
my_prex_or_die "mpirun python metrics_2d.py -i $INPUTDIR -o $METRICS_2D -m $INGVMASK"

my_prex_or_die "python Kz_anom_synthesis.py -i $METRICS_2D -m $INGVMASK -o $PUBDIR"
my_prex_or_die "python metrics_2d_percentiles.py -i $METRICS_2D -m $INGVMASK -o $PERCENTILES_DIR"
my_prex_or_die "python metrics_2d_percentiles_plot.py -i $PERCENTILES_DIR -o $PUBDIR -n $RUNNAME"

# my_prex "python metrics_2d_percentiles_plot_multirun.py" # to plot in https://medeaf.ogs.it/internal-validation/Benchmark/eas7_V10C/multirun_vs_eas8/

OUTPUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/FORCINGS/DELTAT
my_prex_or_die "mpirun python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK"

