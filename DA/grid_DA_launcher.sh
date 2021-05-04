#! /bin/bash

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea


. ./profile.inc

OUTPUTDIR=/galileo/home/userexternal/gcoidess/CODICE_ANNA/preproc/DA/
NAME=DA_grid_5000
MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS24_NRTV7C/meshmask.nc
CASE=2
BASIN=subs
GIB=-4.9

my_prex_or_die "python create_GRID_DA.py -o $OUTPUTDIR -n $NAME -m $MASKFILE -c $CASE -b $BASIN -g $GIB" #-v $VARIABLE -c $COASTNESS -b $BASIN -s $STATISTIC "

