#! /bin/bash

#SBATCH --job-name=valid
#SBATCH -N1 -n 1
#SBATCH --time=2:00:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_all_serial

cd $SLURM_SUBMIT_DIR

. ../profile.inc

REMOTEDIR=/medsea_oce/nrt/dev/medfs_simu_clim/2019/
LOCALDIR=/g100_scratch/userexternal/gbolzon0/V11C/EFAS
KEY=/g100/home/userexternal/gbolzon0/.ssh/g100_cmcc-key

my_prex  " echo \" get -R $REMOTEDIR $LOCALDIR \" | sftp -i $KEY  -P 20022 ogs@dtn01.cmcc.it"
