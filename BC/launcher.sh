#! /bin/bash

#SBATCH --job-name=bc
##SBATCH -N1
#SBATCH -N2
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
##SBATCH --mem=300gb
#SBATCH --mem=30gb
#SBATCH --account=OGS_test2528_0
##SBATCH --partition=lrd_all_serial
#SBATCH --partition=dcgp_usr_prod
##SBATCH --partition=g100_usr_prod

cd $SLURM_SUBMIT_DIR

#source /g100_work/OGS23_PRACE_IT/COPERNICUS/sequence.sh
source /leonardo_work/OGS23_PRACE_IT_0/COPERNICUS/sequence.sh
#source /leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/PYENVS/sequence.sh
export PYTHONPATH=/leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/git/bit.sea/src:$PYTHONPATH


unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false
source ../profile.inc

#MASKFILE=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc
MASKFILE=/leonardo_scratch/fast/RIVRD_OGS/OGSTM-BFM-SETUP/MASKS/med24v1_domain_cfg.nc
MASKFILE_CMCC=/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc
INPUTBC=/leonardo_scratch/fast/RIVRD_OGS/V13C_QUID/INPUTBC/

INPUTDIR=$PWD/INPUT_BC_PROFILES
mkdir -p $INPUTDIR
for var in N1p N3n N5s O2o O3c O3h; do
    cp $INPUTBC/${var}.nc $INPUTDIR
done

# NOTE: 
# O3c_O3h_profiles.py DOES NOT WORK.

#my_prex_or_die "python O3c_O3h_profiles.py -m $MASKFILE -M $MASKFILE_CMCC -o $INPUTDIR"

my_prex_or_die "python write_VPnutrients.py -i $INPUTDIR -m $MASKFILE"
my_prex_or_die "python main.py"

RESTARTSDIR=/leonardo_scratch/fast/RIVRD_OGS/OGSTM-BFM-SETUP/RST/
OUTDIR=/leonardo_scratch/fast/RIVRD_OGS/OGSTM-BFM-SETUP/BC/

cd atlantic
my_prex_or_die "python atlantic_generator.py -i $INPUTBC -o $OUTDIR -m $MASKFILE --rst $RESTARTSDIR"


#cd ../EFAS
#CMCC_MASK=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc
#FORCINGS=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/FORCINGS/READY_FOR_MODEL/
#my_prex_or_die "mpirun python river_generation.py -i $INPUTDIR -o $OUTDIR -m $MASKFILE -M $CMCC_MASK"

