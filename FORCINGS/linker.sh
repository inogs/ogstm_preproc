#! /bin/bash

DIR=/gpfs/scratch/userexternal/gbolzon0/BI_HOURLY/DAILY/ORIG
for I in `ls $DIR/*nc `; do
   filename=`basename $I`
   var=${filename:${#filename}-4:1}
   date8=${filename:${#filename}-18:8}
   ln -fs $I ${var}${date8}-12:00:00.nc
done

ln -fs T20190101-12:00:00.nc T20190101-00:00:00.nc
ln -fs U20190101-12:00:00.nc U20190101-00:00:00.nc
ln -fs V20190101-12:00:00.nc V20190101-00:00:00.nc
ln -fs W20190101-12:00:00.nc W20190101-00:00:00.nc


~/lnk_date 20191231-12:00:00 20200101-12:00:00

