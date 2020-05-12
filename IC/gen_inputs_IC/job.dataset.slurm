#!/bin/bash

#SBATCH --job-name=dtset
#SBATCH -N1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --account=IscrC_REBIOMED 
#SBATCH --partition=gll_all_serial
#SBATCH --error dataset.err
#SBATCH --output dataset.out

cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/scratch/userexternal/vdibiagi/bit.sea-git

date

WRKDIR=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/IC/preproc/IC/gen_inputs_IC
OUTDIR_parent=$WRKDIR

for varName in N1p, N3n, N5s, O2o; do
# period 1
python $WRKDIR/statsProfilesFromDataset.py -p $OUTDIR_parent/ -v $varName -s 1997 -e 2007
done

for varName in N1p, N3n, N5s, O2o, DIC, ALK; do 
#period 2
python $WRKDIR/statsProfilesFromDataset.py -p $OUTDIR_parent/ -v $varName -s 1997 -e 2015
done

date
exit 0



