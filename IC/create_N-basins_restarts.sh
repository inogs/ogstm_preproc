#!/bin/bash
source /g100_work/OGS23_PRACE_IT/ggalli00/VENVZ/sequence_gg.sh

echo computing mean subbasin R3l
python create_R3l_nudging.py

outdir='/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/FILES/'
rstpfx='RST.19950101-00:00:00'
ndgpfx='R3l_yyyy0630-00:00:00.nc'
maskfile='/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022_v6/wrkdir/MODEL/meshmask.nc'
csvfile='/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/SCRIPTS/mean_R3l_from_RRS412.csv'
python create_N-basins_restarts.py -o $outdir -r $rstpfx -n $ndgpfx -m $maskfile -c $csvfile

