import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Step 1 for the generation of a climatology

    Generates a .npy file of the map occurrency, a 2d map of integers, containing
    the number of measurements found in all the daily files for each i,j indexes
    That file will be used to generate the Sat maskfile
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '''ORIG sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG'''
                                )
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                help = '''CHL or KD490, the name of the NetCDF variable in ORIG files'''
                                )
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = True,
                                help = '''CHL_map_occurency.npy'''
                                )
    return parser.parse_args()

args = argument()


from bitsea.commons.Timelist import TimeList
import bitsea.Sat.SatManager as Sat
import numpy as np
from bitsea.commons.utils import addsep

INPUTDIR=addsep(args.inputdir)


TL = TimeList.fromfilenames(None, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')

COUNT = np.zeros((1580,3308),np.int32)
for ifile, filename in enumerate(TL.filelist):
    print "Working on file number ", ifile, " of ", len(TL.filelist)
    fillvalue=-999
    A=Sat.readfromfile(filename,args.var)
    ii=np.isnan(A)
    A[ii] = fillvalue
    filled = A==fillvalue
    goods = ~filled 
    COUNT[goods] = COUNT[goods]+1
    


np.save(args.outfile,COUNT)