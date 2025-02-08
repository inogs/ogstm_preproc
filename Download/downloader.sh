#! /bin/bash

#SBATCH --job-name=MDS
#SBATCH -N1 -n 1
#SBATCH --time=4:00:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_all_serial


cd $SLURM_SUBMIT_DIR

. ../profile.inc

ORIGDIR=
mkdir -p $ORIGDIR

DATASET_LIST="
cmems_obs-oc_med_bgc-plankton_my_l3-multi-1km_P1D
cmems_obs-oc_med_bgc-transp_my_l3-multi-1km_P1D
"
#cmems_obs-oc_med_bgc-plankton_my_l3-olci-300m_P1D
#cmems_obs-oc_med_bgc-transp_my_l3-olci-300m_P1D
#cmems_obs-oc_med_bgc-reflectance_my_l3-multi-1km_P1D"

COMMONS=" -nd -s files --force-download --show-outputnames --overwrite --disable-progress-bar "
for year in $(seq 1999 2023 ) ; do
    for dataset in $DATASET_LIST; do
        my_prex "copernicusmarine get $COMMONS -o ${ORIGDIR}  -i ${dataset}  --filter *${year}/*"
    done
done


