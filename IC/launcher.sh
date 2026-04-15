#! /bin/bash

ORIGMASK=/g100_work/OGS_test2528/Benchmark/SETUP/PREPROC/MASK/meshmask.nc
MASKFILE=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc
INPUTDIR=/g100_work/OGS_test2528/V12C_QUID/PREPROC/IC/RST
OUTPUTDIR=/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/IC/RST
DATE=20220101

mkdir -p $OUTPUTDIR
python interpolator_same_resolution.py -i $INPUTDIR -o $OUTPUTDIR --origmask $ORIGMASK --newmask $MASKFILE -d $DATE
python R8.py -i $OUTPUTDIR -m $MASKFILE -d $DATE-00:00:00


python ICgenerator_fromInputs.py -i $INPUTDIR -o $OUTPUTDIR -m $MASKFILE -d $DATE



