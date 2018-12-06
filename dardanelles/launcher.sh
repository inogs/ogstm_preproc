#! /bin/bash

python u_slice -o dard_slice.npy
python subbasin_analysis.py -o AEG_integrals_24.npy
python monthly_statistics.py -i dard_slice.npy -o figure.png # only statistics

python preserving_ludwig.py -u dard_slice.npy -i AEG_integrals_24.npy -o dard_boundary_values.npy

OUTPUTDIR=/gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/INPUTS_for_MODEL/2/
python dardanelles_generator.py -i dard_boundary_values.npy -o $OUTPUTDIR
