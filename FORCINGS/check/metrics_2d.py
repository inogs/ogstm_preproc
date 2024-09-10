import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates 2d metrics for every daily forcing


    Provided metrics:
    - mld
    - KE_ratio, integrated on 0-200m
    - KE_total, integrated on 0-200m
    - Stratification index defined as integral on 0-200m of Brunt Vaissala frequency * z
    - Vertical eddy diffusivity integrated on 0-200m
    - Vertical eddy diffusivity integrated on 200-500m
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

from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
import numpy as np
from bitsea.commons.Timelist import TimeList
from bitsea.commons.utils import addsep
import seawater as sw
from bitsea.surf import surfaces
from bitsea.commons.layer import Layer
from bitsea.layer_integral.mapbuilder import MapBuilder
from bitsea.commons import netcdf4

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

BathyCells = TheMask.bathymetry_in_cells()

Au = np.zeros((jpk,jpj,jpi),np.float32)
Av = np.zeros((jpk,jpj,jpi),np.float32)
Aw = np.zeros((jpk,jpj,jpi),np.float32)

for k in range(jpk):
    Aw[k,:] = TheMask.area
for k in range(jpk):
    Au[k,:] = TheMask.e2t*TheMask.e3t[k,:]
    Av[k,:] = TheMask.e1t*TheMask.e3t[k,:]


izlev=TheMask.getDepthIndex(200) + 5
iz200=TheMask.getDepthIndex(200) + 2
M  =np.zeros((jpk,jpj,jpi),np.float32)*np.nan
BVF=np.zeros((jpk,jpj,jpi),np.float32)
P  =np.zeros((izlev,jpj,jpi),np.float32)
mask = TheMask.mask_at_level(0)

for ji in range(jpi):
    for jj in range(jpj):
        M[:,jj,ji]=TheMask.zlevels
        P[:,jj,ji]=TheMask.zlevels[:izlev]


layer200 = Layer(0,200)
layer500 = Layer(200,500)

jk_50 = TheMask.getDepthIndex(50)
jk_100 = TheMask.getDepthIndex(100)
jk_150 = TheMask.getDepthIndex(150)

TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", prefix="U")
eps=1.e-08

nFrames=TL.nTimes


FRAMES=range(nFrames)

for iframe in FRAMES[rank::nranks]:
    timestr = TL.Timelist[iframe].strftime("%Y%m%d-%H:%M:%S")
    print(timestr,flush=True)
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

    Eddy_diff_100 = MapBuilder.get_layer_average(K, Layer(100,150))
    Eddy_diff_400 = MapBuilder.get_layer_average(K, Layer(400,500))

    Ratio_Eddy_diff = Eddy_diff_400/Eddy_diff_100

    K_clean = DataExtractor(TheMask, rawdata=K.values ) #copy
    Count_anomalous_Eddy_diff = np.zeros((jpj,jpi),np.float32)
    Count_anomalous_Eddy_diff[~TheMask.mask_at_level(0)] = 1.e+20
    for i in range(jpi):
        for j in range(jpj):
            counter=0
            nCells = BathyCells[j,i]
            if nCells > jk_150: # index of 150m
                Kmin = K.values[jk_50:jk_150,j,i].min()
                if Kmin < 1.e-5:
                    for k in range(jk_150,nCells):
                        if K.values[k,j,i] > 1.e-4:
                            counter+=1
                            K_clean.values[k,j,i] = 1.e-7
                    Count_anomalous_Eddy_diff[j,i] = counter


    Eddy_diff_150 = MapBuilder.get_layer_average(K_clean, Layer(0,150))
    Eddy_diff_500 = MapBuilder.get_layer_average(K_clean, Layer(150,500))


    U[U==0]=eps
    V[V==0]=eps
    W[W==0]=eps
    U[U>1.e+19]=eps
    V[V>1.e+19]=eps
    W[W>1.e+19]=eps
    KE_ratio = W**2*Aw/(U**2*Au + V**2*Av)

    ii=(U==eps) & (V==eps) # taking in account umask and vmask False
    KE_ratio[ii] = 0

    KE_ratio_2d = KE_ratio[:iz200,:,:].sum(axis=0)
    lmask = (KE_ratio_2d>0) & mask
    KE_ratio_2d[lmask] = np.log10(KE_ratio_2d[lmask])
    KE_ratio_2d[~lmask] = 1.e+20


    KE_total = 0.5*(U**2*Au + V**2*Av + W**2*Aw)
    KE_total_2d = KE_total[:iz200,:,:].sum(axis=0)

    KE_vert = 0.5*(W**2*Aw)
    KE_vert_2d = KE_vert[:iz200,:,:].sum(axis=0)

    BVF[:izlev-1,:], _, p_ave = sw.bfrq(S[:izlev,:],T[:izlev,:],P)
    BVF[BVF<0]=0 # because sw.bfrw returns a negative value when S=0, T=0 is encountered as first land point
    De = DataExtractor(TheMask,rawdata=BVF*M)
    Stratification_index = MapBuilder.get_layer_integral(De, layer200)

    Stratification_index[~mask] = 1.e+20
    KE_total_2d[~mask] = 1.e+20
    KE_vert_2d[~mask] = 1.e+20
    mld[~mask] = 1.e+20
    Eddy_diff_150[~mask] = 1.e+20
    Eddy_diff_500[~mask] = 1.e+20

    outfile=OUTPUTDIR + "metrics." + timestr + ".nc"
    netcdf4.write_2d_file(KE_ratio_2d         ,'KE_ratio', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(KE_total_2d         ,'KE_total', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(KE_vert_2d         ,'KE_vert', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Stratification_index,'stratification_index', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(mld                 ,'mld' , outfile,TheMask, compression=True)
    netcdf4.write_2d_file(Eddy_diff_150       ,'Ved_150', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Eddy_diff_500       ,'Ved_500', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Ratio_Eddy_diff     ,'Ved_ratio', outfile, TheMask, compression=True)
    netcdf4.write_2d_file(Count_anomalous_Eddy_diff,'anom_counter', outfile, TheMask, compression=True)

