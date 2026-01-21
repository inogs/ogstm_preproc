#!/bin/bash
#SBATCH -N1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --mem=256gb
#SBATCH -o dr.%j.out
#SBATCH -e dr.%j.err
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=OGS_test2528
##SBATCH --account=OGS23_Tbit10

cd $SLURM_SUBMIT_DIR

source /leonardo_work/OGS23_PRACE_IT_0/ggalli00/NECCTON/SCRIPTS/VENV/sequence.sh

echo START: `date`

echo degrading restarts

gitdir=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/git/ogstm_preproc/FORCINGS/degradation
cd $gitdir
git checkout 15-forcing-degradation

mpirun -n 8 python degrade_restart.py -y degrade_restart.yaml

echo END: `date`
