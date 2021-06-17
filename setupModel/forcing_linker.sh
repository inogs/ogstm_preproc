#! /bin/bash

# ####################################################################
# ogstm uses forcing names like T20160101-00:00:00.nc 
# ln -s  assw_v4_1d_20151231_20160101_grid_T.nc T20160101-00:00:00.nc
######################################################################


DIR=/gpfs/scratch/userexternal/gcoidess/TRANSITION_PREOPERATIONAL
for I in `ls $DIR/mfs[12]*nc `; do
   filename=`basename $I`
   var=${filename:${#filename}-4:1}
   date8=${filename:${#filename}-18:8}
   ln -fs $I ${var}${date8}-12:00:00.nc
done

ln -fs T20190101-12:00:00.nc T20190101-00:00:00.nc
ln -fs U20190101-12:00:00.nc U20190101-00:00:00.nc
ln -fs V20190101-12:00:00.nc V20190101-00:00:00.nc
ln -fs W20190101-12:00:00.nc W20190101-00:00:00.nc

