#! /bin/bash

MASKFILE=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc
INPUTBC=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/BC/inputs

python write_VPnutrients.py -i $INPUTBC -m $MASKFILE
python main.py

RESTARTSDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/IC/RST_2018
OUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/BC/out

cd atlantic
python atlantic_generator.py -i $INPUTBC -o $OUTDIR -m $MASKFILE --rst $RESTARTSDIR

cd ..
CMCC_MASK=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask_CMCC.nc
FORCINGS=/g100_scratch/userexternal/gbolzon0/V9C/2019/FORCINGS/
python Po/online/po_generation -i $FORCINGS -o $OUTDIR -M $CMCC_MASK -m $MASKFILE

