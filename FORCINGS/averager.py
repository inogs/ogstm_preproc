import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Executes ncea in parallel of forcing files, U, V, W T
    Works on files like U20190101-12:00:00.nc
    Performs time average in order to get files with low frequency (4h, 6h)
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'Directory with High Frequency files')

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = 'Directory with Low Frequency files')
    parser.add_argument(   '--frequency', '-f',
                                type = str,
                                required = True,
                                help = 'Output frequency in h')
        
    return parser.parse_args()

args = argument()


from commons.Timelist import TimeList
from commons import timerequestors
from commons import genUserDateList as DL
import os
from commons.utils import addsep
from datetime import datetime, timedelta

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1

INPUTDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outputdir)

dateformat="%Y%m%d-%H:%M:%S"

#Start_time="20210101-00:00:00"
#End_time = "20220101-00:00:00"
output_frequency_h= int(args.frequency)

#start=datetime.strptime(Start_time,dateformat) + timedelta(hours=output_frequency_h/2)

#OUT_TIMES = DL.getTimeList(start.strftime(dateformat),End_time, hours=output_frequency_h)



for var in ["U", "V", "W", "T"]:
    TL = TimeList.fromfilenames(None, INPUTDIR, var + "*nc",prefix=var)
    d0=TL.Timelist[0]
    dm1=TL.Timelist[-1]
    End_time=dm1.strftime(dateformat)
    start=datetime(d0.year,d0.month,d0.day) + timedelta(hours=output_frequency_h/2)
    OUT_TIMES = DL.getTimeList(start.strftime(dateformat),End_time, hours=output_frequency_h)


    for d in OUT_TIMES[rank::nranks]:
        req = timerequestors.Hourly_req(d.year,d.month,d.day, d.hour, delta_hours=output_frequency_h)
        ii,w = TL.select(req)
        files_string = ' '.join([str(k) for k in TL.filelist])
        outfile = OUTDIR + var + d.strftime(dateformat) + ".nc"
        print("rank", rank, outfile, flush=True)
        command = "ncea %s -O %s" %( files_string, outfile)
        os.system(command)
