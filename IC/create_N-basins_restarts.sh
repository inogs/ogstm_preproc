#!/bin/bash

. ../profile.inc
source /g100_work/OGS23_PRACE_IT/ggalli00/VENVZ/sequence_gg.sh

csvfile=$PWD/mean_R3l_from_RRS412.csv
my_prex_or_die "python create_R3l_nudging.py -o $csvfile"

outdir=$PWD/
rstpfx=RST.19950101-00:00:00
ndgpfx=R3l_yyyy0630-00:00:00.nc
maskfile=/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc

my_prex_or_die "python create_N-basins_restarts.py -o $outdir -r $rstpfx -n $ndgpfx -m $maskfile -c $csvfile"

