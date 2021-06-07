#! /bin/bash

# ####################################################################
# ogstm uses forcing names like T20160101-00:00:00.nc 
# ln -s  assw_v4_1d_20151231_20160101_grid_T.nc T20160101-00:00:00.nc
######################################################################

DIR=/gpfs/scratch/userexternal/plazzari/eas_v12/FORCINGS/UNZIPPED
for I in `ls $DIR ` ; do
  ogstm_name=${I:35:1}${I:21:8}-00:00:00.nc
  ln -s ${DIR}/${I} $ogstm_name
done

