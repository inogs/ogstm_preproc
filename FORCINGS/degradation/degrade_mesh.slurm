#!/bin/bash
#SBATCH -N1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem=32gb
#SBATCH -o dm.%j.out
#SBATCH -e dm.%j.err
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=OGS_test2528
##SBATCH --account=OGS23_Tbit10

cd $SLURM_SUBMIT_DIR

source /leonardo_work/OGS23_PRACE_IT_0/ggalli00/NECCTON/SCRIPTS/VENV/sequence.sh

echo DEGRADING MESH
echo START: `date`

gitdir=/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/REGRID/git/ogstm_preproc/FORCINGS/degradation
cd $gitdir
git checkout 15-forcing-degradation

python mesh_degradation.py -y degrade_mesh.yaml

echo END: `date`
