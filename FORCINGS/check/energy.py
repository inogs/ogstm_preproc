import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''

    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/gpfs/work/OGS_prod_0/OPA/V5C/devel/wrkdir/2/MODEL/FORCINGS/')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc')

    return parser.parse_args()

args = argument()

from commons.mask import Mask
from commons.dataextractor import DataExtractor
import numpy as np
from commons.Timelist import TimeList
from commons.utils import addsep
import seawater as sw
from surf import surfaces
from commons.layer import Layer
from layer_integral.mapbuilder import MapBuilder
from commons import netcdf4

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outdir)
TheMask=Mask(args.maskfile)

jpk, jpj, jpi = TheMask.shape
izlev=TheMask.getDepthIndex(200) + 5 #350
M=np.zeros((jpk,jpj,jpi),np.float32)*np.nan
BVF=np.zeros((jpk,jpj,jpi),np.float32)
mask = TheMask.mask_at_level(0)

for ii in range(jpi):
    for jj in range(jpj):
        M[:,jj,ii]=TheMask.zlevels
P=np.transpose([[TheMask.zlevels[:izlev]]* jpj] * jpi )

layer200 = Layer(0,200)
layer500 = Layer(200,500)


TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", prefix="U")
eps=1.e-08

nFrames=TL.nTimes


FRAMES=range(nFrames)

for iframe in FRAMES[rank::nranks]:
    timestr = TL.Timelist[iframe].strftime("%Y%m%d-%H:%M:%S")
    print(timestr)
    filenameU=INPUTDIR + "U" + timestr + ".nc"
    filenameV=INPUTDIR + "V" + timestr + ".nc"
    filenameW=INPUTDIR + "W" + timestr + ".nc"
    filenameT=INPUTDIR + "T" + timestr + ".nc"

    U=DataExtractor(TheMask,filenameU,"vozocrtx").values
    V=DataExtractor(TheMask,filenameV,"vomecrty").values
    W=DataExtractor(TheMask,filenameW,"vovecrtz").values
    K=DataExtractor(TheMask,filenameW,"votkeavt")
    T=DataExtractor(TheMask,filenameT,"votemper").values
    S=DataExtractor(TheMask,filenameT,"vosaline").values
    mld=np.abs(DataExtractor(TheMask,filenameT,"somxl010").values)
    
    MLD = surfaces.mld(T, TheMask)
    Eddy_diff_200 = MapBuilder.get_layer_average(K, layer200)
    Eddy_diff_500 = MapBuilder.get_layer_average(K, layer500)


    U[U==0]=eps
    V[V==0]=eps
    W[W==0]=eps
    U[U>1.e+19]=eps
    V[V>1.e+19]=eps
    W[W>1.e+19]=eps
    KE_ratio = W**2/(U**2 + V**2)

    ii=(U==eps) & (V==eps) # taking in account umask and vmask False
    KE_ratio[ii] = 0

    De = DataExtractor(TheMask,rawdata=KE_ratio)
    KE_ratio_2d = MapBuilder.get_layer_integral(De,layer200)
    lmask = (KE_ratio_2d>0) & mask
    KE_ratio_2d[lmask] = np.log10(KE_ratio_2d[lmask])
    KE_ratio_2d[~lmask] = 1.e+20


    KE_total = 0.5*(U**2 + V**2 + W**2)
    De = DataExtractor(TheMask,rawdata=KE_total)
    KE_total_2d = MapBuilder.get_layer_integral(De,layer200)

    BVF[:izlev-1,:], _, p_ave = sw.bfrq(S[:izlev,:],T[:izlev,:],P)
    BVF[BVF<0]=0 # because sw.bfrw returns a negative value when S=0, T=0 is encountered as first land point
    De = DataExtractor(TheMask,rawdata=BVF*M)
    Stratification_index = MapBuilder.get_layer_integral(De, layer200)

    Stratification_index[~mask] = 1.e+20
    KE_total_2d[~mask] = 1.e+20
    mld[~mask] = 1.e+20
    Eddy_diff_200[~mask] = 1.e+20
    Eddy_diff_500[~mask] = 1.e+20

    outfile=OUTPUTDIR + "metrics." + timestr + ".nc"
    netcdf4.write_2d_file(KE_ratio_2d         ,'KE_ratio', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(KE_total_2d         ,'KE_total', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Stratification_index,'stratification_index', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(mld                 ,'mld' , outfile,TheMask, compression=True)
    netcdf4.write_2d_file(Eddy_diff_200       ,'Evd_200', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Eddy_diff_500       ,'Evd_500', outfile, TheMask, compression=True)
