#! /bin/bash

#SBATCH --job-name=WGET
#SBATCH -N1 -n 1
#SBATCH --time=4:00:00
#SBATCH --account=OGS_dev
#SBATCH --partition=gll_all_serial

cd $SLURM_SUBMIT_DIR

. ../profile.inc

PRODUCT=OCEANCOLOUR_MED_BGC_L3_MY_009_143

#cmems_obs-oc_med_bgc-optics_my_l3-multi-1km_P1D
DATASET_LIST="
cmems_obs-oc_med_bgc-reflectance_my_l3-olci-300m_P1D
cmems_obs-oc_med_bgc-plankton_my_l3-multi-1km_P1D
cmems_obs-oc_med_bgc-transp_my_l3-multi-1km_P1D
cmems_obs-oc_med_bgc-plankton_my_l3-olci-300m_P1D
cmems_obs-oc_med_bgc-transp_my_l3-olci-300m_P1D
cmems_obs-oc_med_bgc-reflectance_my_l3-multi-1km_P1D"

for dataset in $DATASET_LIST; do 

   my_prex "wget --ftp-user=MED_OGS_TRIESTE_IT --ftp-password=NEdifupa -m ftp://my.cmems-du.eu/Core/${PRODUCT}/${dataset}"

done
