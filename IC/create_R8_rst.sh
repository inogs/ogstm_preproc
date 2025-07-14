#!/bin/bash

# also adds R8 to the ATL BC file

module load nco

INDIR_RS=/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022/wrkdir/IC
OUTDIR_RS=/g100_work/OGS23_PRACE_IT/ggalli00/OGSTM-BFM/OGSTM-BFM_setup/RESTARTS

INDIR_BC=/g100_work/OGS_devC/Benchmark/HC/wrkdir/MODEL/BC
OUTDIR_BC=/g100_work/OGS23_PRACE_IT/ggalli00/OGSTM-BFM/OGSTM-BFM_setup/BCs

constituents=(c n p s)

#rsync -aP $INDIR_BC/ATL_yyyy0630-00:00:00.nc $OUTDIR_BC/ATL_yyyy0630-00:00:00.nc

RST_HEAD=RST.19990101-00:00:00.

for currc in "${constituents[@]}"; do

  echo $currc

#  # make restart
  rsync -aP $INDIR_RS/${RST_HEAD}R6${currc}.nc $OUTDIR_RS/${RST_HEAD}R8${currc}.nc
  ncap2 -O -s "TRNR8${currc}=TRNR6${currc}*0.1" $OUTDIR_RS/${RST_HEAD}R8${currc}.nc $OUTDIR_RS/${RST_HEAD}R8${currc}.nc
  ncap2 -O -s "where(TRNR8${currc} > 1.e+10) TRNR8${currc}=1.e+20" $OUTDIR_RS/${RST_HEAD}R8${currc}.nc $OUTDIR_RS/${RST_HEAD}R8${currc}.nc
  ncks -C -O -x -v TRNR6${currc} $OUTDIR_RS/${RST_HEAD}R8${currc}.nc $OUTDIR_RS/${RST_HEAD}R8${currc}.nc

  # make ATL BC
#  ncap2 -O -s "R8${currc}=R6${currc}*0.1" $OUTDIR_BC/ATL_yyyy0630-00:00:00.nc $OUTDIR_BC/ATL_yyyy0630-00:00:00.nc
#  ncap2 -O -s "where(R8${currc} > 1.e+10) R8${currc}=1.e+20" $OUTDIR_BC/ATL_yyyy0630-00:00:00.nc $OUTDIR_BC/ATL_yyyy0630-00:00:00.nc

done

