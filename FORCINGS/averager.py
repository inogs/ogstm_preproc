import argparse
from bitsea.utilities.argparse_types import existing_dir_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Executes ncea in parallel of forcing files, U, V, W T
    Works on files like U20190101-12:00:00.nc
    Performs time average in order to get files with low frequency (4h, 6h)
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = 'Directory with High Frequency files')

    parser.add_argument(   '--outputdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = 'Directory with Low Frequency files')
    parser.add_argument(   '--frequency', '-f',
                                type = str,
                                required = True,
                                help = 'Output frequency in h')
        
    return parser.parse_args()

args = argument()


from bitsea.commons.Timelist import TimeList
from bitsea.commons import timerequestors
from bitsea.commons import genUserDateList as DL
import os
from datetime import datetime, timedelta
import mpi4py
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator
comm = get_mpi_communicator()
rank  = comm.Get_rank()
nranks =comm.size

INPUTDIR = args.inputdir
OUTDIR = args.outputdir

dateformat="%Y%m%d-%H:%M:%S"

output_frequency_h= int(args.frequency)



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
        files_string = ' '.join([str(TL.filelist[k]) for k in ii])
        outfile = OUTDIR / '{}{}.nc'.format(var,d.strftime(dateformat))
        print("rank", rank, outfile, flush=True)
        command = "ncea %s -O %s" %( files_string, outfile)
        os.system(command)
