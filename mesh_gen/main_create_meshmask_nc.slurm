#!/bin/bash

#SBATCH --job-name=mesh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem=24gb
#SBATCH --account=OGS23_PRACE_IT_0
##SBATCH --partition=dcgp_usr_prod
#SBATCH --partition=lrd_all_serial 

cd $SLURM_SUBMIT_DIR

source /leonardo_work/OGS23_PRACE_IT_0/COPERNICUS/sequence.sh

mesh_in=/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc
mesh_out=/leonardo_work/RIVRD_OGS/OGSTM-BFM-SETUP/MASKS/med24v1_domain_cfg.nc

gitdir=/leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/git/ogstm_preproc/
cd $gitdir/mesh_gen
#git checkout rivierade

#python main.py
python create_meshmask_nc.py -i $mesh_in -o $mesh_out
