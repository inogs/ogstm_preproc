import numpy as np
from bitsea.commons import netcdf4
from bitsea.commons.Timelist import TimeList

INPUTDIR="/gpfs/scratch/userexternal/gcoidess/EAS6_2019_PO/"

riverinput="/gpfs/scratch/userexternal/gbolzon0/PO/runoff_1d_nomask_y2019.nc"

S = netcdf4.readfile(riverinput, 's_river')
R = netcdf4.readfile(riverinput, 'sorunoff')


ii = S[0,:,:]==18
J,I = np.nonzero(ii)

nPoints=ii.sum()

TL = TimeList.fromfilenames(None, INPUTDIR, "T*.nc", prefix="assw_drpo_2019_1d_", dateformat="%Y%m%d")

nFrames = TL.nTimes
PO_RUNOFF=np.zeros((nFrames,nPoints), np.float32)
INPUT_RUNOFF=np.zeros((nFrames,nPoints), np.float32)


for iFrame, filename in enumerate(TL.filelist):
    print iFrame
    nemo_Runoff=netcdf4.readfile(filename, "sorunoff")[0,:]
    PO_RUNOFF[iFrame,:] = nemo_Runoff[ii]
    INPUT_RUNOFF[iFrame,:] = R[iFrame,ii]
    
np.save('nemo_runoff.npy',PO_RUNOFF)
np.save('input_runoff.npy',INPUT_RUNOFF)


LIST=['Po Volano',
      'Po Goro',
      'Po Gnocca',
      'Po Bocca Tolle',
      'Po Bastimento',
      'Po Scirocco + Po Bonifazio',
      'Po Dritta',
      'Po Tramontana',
      'Po Maistra ',
      'Po Levante']



import pylab as pl


for k in range(10):
    fig,ax = pl.subplots()
    ax.plot(TL.Timelist,PO_RUNOFF[:,k],'b-',label='ref')
    ax.plot(TL.Timelist,INPUT_RUNOFF[:,k],'r.',label='model')
    ax.legend()
    outfile="P%d.png" %(k)
    fig.savefig(outfile)
    pl.close(fig)


RIVER_CONCENTRATION={'O2o':250  ,  #mmol/m3
                     'N1p':2.572,  #mmol/m3
                     'N3n':150  ,  #mmol/m3
                     'N5s':150  ,  #mmol/m3
                     'O3c':33225,  #  mg/m3
                     'O3h':2800 }  #mmol/m3

import config as conf
from bclib.river import river
conf.file_river = 'Perseus-4.6_39rivers_mesh24.xlsx'

PERSEUS = river(conf)
PERSEUS.modularize(conf)
YEARS = np.arange(2000,2011)
CLIM={}
for sheet in conf.river_data_sheet:
    print sheet
    SUM = 0.
    for year in YEARS:
        year_str = str(year)
        annual_contribution=PERSEUS.xls_data[sheet][year_str][3:13].sum()
        print year, annual_contribution
        SUM += annual_contribution
    CLIM[sheet] = SUM/len(YEARS)
