#!/bin/bash

#SBATCH --job-name=DREScopy
#SBATCH -N1 -n 1
#SBATCH --time=4:00:00
#SBATCH --account=OGS_dev_1
#SBATCH --partition=gll_all_serial

cd $SLURM_SUBMIT_DIR

. ../profile.inc
echo "Archiving on DRES"

INPUTDIR=/gpfs/scratch/userexternal/gcoidess/TRANSITION/SIM_PHYS_daily
OUTDIR=/gss/gss_work/DRES_OGS_BiGe/ateruzzi/EAS/EIS_V8C/SIM_PHYS_daily

cd $INPUTDIR 
for var in T U V W ; do
   for m in `seq 1 12 `; do
      month=`printf %02d $m`
      tarfile=$OUTPUTDIR/${var}2019${month}.tar
      my_prex_or_die  "tar -cf $tarfile simw_eas6_2ht_1d_2019${month}??_grid_${var}.nc "
   done
done

