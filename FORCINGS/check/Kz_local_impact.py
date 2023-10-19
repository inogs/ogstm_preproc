import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
from commons.mask import Mask
from commons.dataextractor import DataExtractor

TheMask=Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc')
DIR1="/g100_scratch/userexternal/gbolzon0/V10C/TEST_FORCINGS_24h/wrkdir/MODEL/AVE_FREQ_1/"
DIR2="/g100_scratch/userexternal/gbolzon0/V10C/run1.1/wrkdir/MODEL/AVE_FREQ_1/"
OUTDIR=""


def plot_profile(i,j, outfile):
    k = B[j,i]
    fig,ax=pl.subplots()
    fig.set_size_inches(5,10)
    
    z=TheMask.zlevels[:k]
    ax.plot(V_orig[:k,j,i], z, 'b', label='orig')
    ax.plot(V_corr[:k,j,i], z, 'r', label='corr')
    ax.grid(True)
    ax.invert_yaxis()
    ax.legend()
    ax.set_ylabel('z (m)')
    fig.suptitle("%s  i=%d, j=%d" %(var,j,i))
    fig.savefig(outfile)


date="20190413"

I = [410,904,370,663]
J = [181,129,195,214]

i,j = 410, 181
i,j = 904, 129
i,j = 370, 195
i,j = 663, 214 # la correzione non ha preso in z<150
for var in ["N1p", "N3n", "O2o"]:

    outfile="%s/worstcase_%s.%s.png" %(OUTDIR,date,var)

    filename="ave." + date + "-12:00:00." + var + ".nc"
    V_orig=DataExtractor(TheMask, DIR1+filename, var ).values
    V_corr=DataExtractor(TheMask, DIR2+filename, var ).values
    B = TheMask.bathymetry_in_cells()

    D=np.abs(V_corr-V_orig)
    K,J,I=np.where(D==D.max())
    k,j,i = K[0],J[0],I[0]

    plot_profile(i, j, outfile)
    for k in range(len(I)):
        i=I[k]
        j=J[k]
        outfile="%s/%d_%d_%s.%s.png" %(OUTDIR,i,j,date,var)
        plot_profile(i, j, outfile)
    # k = B[j,i]
    # fig,ax=pl.subplots()
    # fig.set_size_inches(5,10)
    #
    # z=TheMask.zlevels[:k]
    # ax.plot(V_orig[:k,j,i], z, 'b', label='orig')
    # ax.plot(V_corr[:k,j,i], z, 'r', label='corr')
    # ax.grid(True)
    # ax.invert_yaxis()
    # ax.legend()
    # ax.set_ylabel('z (m)')
    # fig.suptitle("%s  i=%d, j=%d" %(var,j,i))
    # fig.savefig(outfile)
