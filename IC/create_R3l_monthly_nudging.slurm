#!/bin/bash
#SBATCH --job-name=r3l
#SBATCH -N 1 -n 1
#SBATCH --time=01:00:00
#SBATCH --mem=100gb
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod

cd $SLURM_SUBMIT_DIR

source /g100_work/OGS23_PRACE_IT/ggalli00/VENVZ/sequence_gg.sh

python create_R3l_monthly_nudging.py
