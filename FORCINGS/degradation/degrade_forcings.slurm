#!/bin/bash
#SBATCH -N2
#SBATCH --ntasks-per-node=8
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

echo degrading phys forcings

gitdir=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/git/ogstm_preproc/FORCINGS/degradation
cd $gitdir
git checkout 15-forcing-degradation

#yamlfile=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/degrade_forcings.yaml

# must have ran this beforehand
#python degrade_mesh.py

mpirun -n 16 python degrade_forcings.py -y degrade_forcings.yaml
#mpirun -n 16 python degrade_forcings.py -y $yamlfile
#mpirun -n 16 python degrade_forcings.py -y $yamlfile

echo END: `date`
