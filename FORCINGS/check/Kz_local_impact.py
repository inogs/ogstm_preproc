import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
from commons.mask import Mask
from commons.dataextractor import DataExtractor
import os
from commons import netcdf4
from datetime import datetime, timedelta

TheMask=Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc')
DIR1="/g100_scratch/userexternal/gbolzon0/V10C/run1.2/wrkdir/MODEL/AVE_FREQ_1/"
DIR2="/g100_scratch/userexternal/gbolzon0/V10C/run1.3/wrkdir/MODEL/AVE_FREQ_1/"
FORCINGSDIR="/g100_scratch/userexternal/gbolzon0/V10C/run1.2/wrkdir/MODEL/FORCINGS/"

OUTDIR="/g100_work/OGS_devC/Benchmark/pub/gbolzon/eas8/24h/correction_effect"


def plot_profile(i,j, outfile):
    k = B[j,i]
    fig,ax=pl.subplots()
    fig.set_size_inches(5,10)
    
    z=TheMask.zlevels[:k]
    ax.plot(V_orig[:k,j,i], z, 'b', label='orig')
    ax.plot(V_corr[:k,j,i], z, 'r', label='corr')
    ax.invert_yaxis()
    ax.grid(True)
    ylim = ax.get_ylim()
    yTicks=ax.get_yticks().tolist()
    yTicks.append(150)
    ax.set_yticks(np.sort(yTicks))

    ax.legend()
    ax.set_ylabel('z (m)')
    fig.suptitle("%s  i=%d, j=%d" %(var,i,j))
    ax.set_ylim(ylim)
    fig.savefig(outfile)
    pl.close(fig)


date="20190215"
os.system('mkdir -p %s/%s' %(OUTDIR,date))

I_survey_points = [410,904,370,663]
J_survey_points = [181,129,195,214]

I_survey_points = [160,336,856,663]
J_survey_points = [137,187,54,185]

pl.close('all')

for var in ["N1p", "N3n", "O2o"]:

    outfile="%s/%s/worstcase_%s.%s.png" %(OUTDIR,date,date,var)
    print(outfile)

    filename="ave." + date + "-12:00:00." + var + ".nc"
    V_orig=DataExtractor(TheMask, DIR1+filename, var ).values
    V_corr=DataExtractor(TheMask, DIR2+filename, var ).values
    B = TheMask.bathymetry_in_cells()

    D=np.abs(V_corr-V_orig)
    K,J,I=np.where(D==D.max())
    k,j,i = K[0],J[0],I[0]
    k = B[j,i]

    plot_profile(i, j, outfile)
 
    dateprev=(datetime.strptime(date,'%Y%m%d')-timedelta(days=1)).strftime("%Y%m%d")
    for d in [date, dateprev]:
        outfileKz="%s/%s/Kz_worst.%s_%s.png" %(OUTDIR,date,var, d)
        print(outfileKz)
        forcings_w="%sW%s-12:00:00.nc" %(FORCINGSDIR,d)
        Kz = netcdf4.readfile(forcings_w, 'votkeavt')[0,:]
        fig,ax=pl.subplots()
        fig.set_size_inches(5,10)
        z=TheMask.zlevels[:k]
        ax.plot(Kz[:k,j,i+222], z, 'b', label='VED')
        ax.invert_yaxis()
        ylim = ax.get_ylim()
        yTicks=ax.get_yticks().tolist()
        yTicks.append(150)
        ax.set_yticks(np.sort(yTicks))
        ax.grid(True)
        ax.legend()
        ax.set_ylabel('Kz (m2/s)')
        fig.suptitle("%s  i=%d, j=%d" %('Kz',i,j))
        ax.set_ylim(ylim)
        fig.savefig(outfileKz)


    for k in range(0): #len(I_survey_points)):
        i=I_survey_points[k]
        j=J_survey_points[k]
        outfile="%s/%s/%d_%d_%s.%s.png" %(OUTDIR,date,i,j,date,var)
        print(outfile)
        plot_profile(i, j, outfile)

