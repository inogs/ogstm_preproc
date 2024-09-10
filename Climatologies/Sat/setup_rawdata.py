import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Step 3 of a generation of a climatology.
    Generates 25 raw datafiles with all the timeseries of the Sat dataset
    Each file refers to a subdomain. Here the domain decomposition is hardcoded 5x5
    The information content of raw data files is the same of ORIG dir
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '''ORIG sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '''CHL_1km_meshmask.nc'''
                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' path of the out raw data dir'''
                                )
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                help = '''CHL or KD490'''
                                )
    return parser.parse_args()

args = argument()

import netCDF4
import numpy as np
from bitsea.commons.Timelist import TimeList
import bitsea.Sat.SatManager as Sat
import dom_dec
from bitsea.commons.utils import addsep

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
    print "MyId = ", rank
except:
    rank   = 0
    nranks = 1
    isParallel = False



INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)

TL = TimeList.fromfilenames(None, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')
nFrames = TL.nTimes

D=netCDF4.Dataset(args.maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 5
nproc_j = 5
PROCESSES = np.arange(nproc_i*nproc_j)


for ip in PROCESSES[rank::nranks]:
    (jproc, iproc) = divmod(ip,nproc_i)
    outfile = OUTDIR + "%s_raw_data_%d_%d" %(args.var, jproc,iproc)

    I_start,I_end, J_start, J_end = dom_dec.dom_dec(iproc, jproc, jpiglo, jpjglo, nproc_i, nproc_j)
    tmask = tmask_glo[J_start:J_end, I_start:I_end]
    print "rank ", ip, "processes " , outfile, iproc, jproc, I_start,I_end, J_start, J_end

    Nwaterpoints = tmask.sum()
    print "Nwaterpoints", Nwaterpoints, ip
    if (Nwaterpoints ==0) : continue

    RAW_DATA = np.zeros((Nwaterpoints,nFrames),np.float32)
    for ifile, filename in enumerate(TL.filelist):
        print "reading", ip, ifile
        A=Sat.readfromfile(filename,args.var)
        A=A[J_start:J_end, I_start:I_end]
        RAW_DATA[:,ifile] = A[tmask]

    np.save(outfile,RAW_DATA)







