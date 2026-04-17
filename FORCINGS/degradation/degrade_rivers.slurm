#!/bin/bash
#SBATCH -N1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --mem=256gb
#SBATCH -o df.%j.out
#SBATCH -e df.%j.err
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=OGS_test2528
##SBATCH --account=OGS23_Tbit10

# takes about 18 min to process 1y forcings (T, U, V, W) on 16 cpus (2 nodes)
# node needs to be under-populated, else OOM

cd $SLURM_SUBMIT_DIR

source /leonardo_work/OGS23_PRACE_IT_0/ggalli00/NECCTON/SCRIPTS/VENV/sequence.sh

echo START: `date`

echo degrading river forcings

gitdir=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/git/ogstm_preproc/FORCINGS/degradation
cd $gitdir
git checkout main


mpirun -n 32 python degrade_rivers.py -y degrade_rivers.yaml

echo END: `date`
