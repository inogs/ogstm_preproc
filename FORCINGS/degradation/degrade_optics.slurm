#!/bin/bash
#SBATCH -N2
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
#SBATCH --mem=256gb
#SBATCH -o do.%j.out
#SBATCH -e do.%j.err
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=OGS_test2528
##SBATCH --account=OGS23_Tbit10

# takes about 15 min to process 1y (atm files only) on 96 cpus (2 nodes)

cd $SLURM_SUBMIT_DIR

source /leonardo_work/OGS23_PRACE_IT_0/ggalli00/NECCTON/SCRIPTS/VENV/sequence.sh

echo START: `date`

echo degrading phys forcings

gitdir=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/git/ogstm_preproc/FORCINGS/degradation
cd $gitdir
git checkout 15-forcing-degradation

# must have ran this beforehand
#python degrade_mesh.py

#yamlfile=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/degrade_optics.yaml

mpirun -n 96 python degrade_optics.py -y degrade_optics.yaml
#mpirun -n 96 python degrade_optics.py -y $yamlfile

echo END: `date`
