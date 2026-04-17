import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')


    parser.add_argument(   '--variable', '-v',
                                type = str,
                                default = None,
                                required = True,
                                help = '''model variable''')
    return parser.parse_args()


args = argument()

from bitsea.commons.utils import addsep
import numpy as np
import pandas as pd
from bitsea.commons import timerequestors
from bitsea.instruments import superfloat
from bitsea.instruments import superfloat as bio_float
from bitsea.instruments.var_conversions import FLOATVARS
from bitsea.basins import V2 as OGS
from bitsea.basins.basin import ComposedBasin
from bitsea.commons.mask import Mask
import sys
import matplotlib.pyplot as plt


SAVE_CSV_AND_PNG=False

def plot_line_profiles(df, z_interp, namesub, varmod='O2o'):
    # Crea la figura e i subplot
    fig, axs = plt.subplots(1, 1, figsize=(10, 6))

    # Primo plot: Valori delle righe (numero di osservazione) sull'asse x, profondità (colonne) sull'asse y
    plt.plot(df.values, z_interp, alpha=0.3)
    plt.gca().invert_yaxis()  # Inverte l'asse y per simulare la profondità crescente verso il basso
    plt.ylabel(varmod)
    plt.xlabel("depth (m)")
    plt.title("All profiles " + namesub)

    # Griglia e sfondo
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
    plt.gca().set_facecolor('silver')

    # Ottimizza la disposizione
    plt.tight_layout()
    plt.savefig('FIGURE/' + namesub+ '_'+ varmod + '.png'  )
    df.to_csv('FIGURE/'   + namesub+ '_'+ varmod + '.csv')

OUTDIR     = addsep(args.outdir)
varmod     = args.variable

# INIT
TheMask=Mask.from_file("/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022/wrkdir/MASKS/meshmask.nc")
z_interp= TheMask.zlevels


if OGS.atl in OGS.Pred.basin_list:
  OGS.Pred.basin_list.remove(OGS.atl) # tolgo Atlantic buffer
else: pass

SUBS    = OGS.Pred.basin_list[:]
CLIM    = np.zeros((len(SUBS),  len(z_interp) ), np.float32)*np.nan
STD     = np.zeros((len(SUBS),  len(z_interp) ), np.float32)*np.nan
print('_________________start__________________')


TI=timerequestors.TimeInterval(
    starttime='20120101',
    endtime='20240101',
    dateformat='%Y%m%d',
)

SUB_COUNT = 0
for ISUB in SUBS:
       print('_____________ '+ str(ISUB)  +' _____________')

       Profilelist     = superfloat.FloatSelector(FLOATVARS[varmod],TI, ISUB )

       if (ISUB.name == "ion1"):
           isub = ComposedBasin('ion4', [OGS.swm2, OGS.ion2, OGS.tyr2], 'Neighbors of ion1')
           Profilelist     = superfloat.FloatSelector(FLOATVARS[varmod],TI, isub )
       elif (ISUB.name == "tyr1"):
           isub = ComposedBasin('supertyr', [OGS.tyr2,OGS.tyr1], 'tyr1and2')
           Profilelist     = superfloat.FloatSelector(FLOATVARS[varmod],TI, isub)
       else:
           Profilelist     = superfloat.FloatSelector(FLOATVARS[varmod],TI, ISUB )
       if Profilelist:
         print("number of profiles used")  
         print(len(Profilelist))  
         SERV_VAR = np.full((len(Profilelist), len(z_interp)), np.nan, dtype=np.float32)
         
         ICONT=0
         for PROFILE in Profilelist:
            Pres, Profile, Qc = PROFILE.read(FLOATVARS[varmod])
            Profile_interp = np.interp(z_interp, Pres, Profile)
            SERV_VAR[ICONT, :] = Profile_interp
            ICONT+=1
       
       if SAVE_CSV_AND_PNG:
          df =pd.DataFrame(SERV_VAR)
          df=df.T
          namesub=ISUB.name
          plot_line_profiles(df, z_interp, namesub )


       serv_P =  np.nanmean(SERV_VAR,axis=0) # allprofiles | aLB [0] | ALL Z [:]
       serv_S =  np.nanstd(SERV_VAR,axis=0)
       CLIM[SUB_COUNT, :] = serv_P
       STD[SUB_COUNT, :]  = serv_S
       SUB_COUNT+=1

#mean
import netCDF4 
outfile='yr_Avg_'+varmod+'.nc'
ncOUT = netCDF4.Dataset(outfile,'w')
ncOUT.createDimension('nsub',len(SUBS) )
ncOUT.createDimension('nav_lev',  len(z_interp) )
ncvar=ncOUT.createVariable(varmod,'f',('nsub','nav_lev'))
ncvar[:]=(CLIM)
ncOUT.close()

#std
outfile='yr_Std_'+varmod+'.nc'
ncOUT = netCDF4.Dataset( outfile,'w')
ncOUT.createDimension('nsub',len(SUBS) )
ncOUT.createDimension('nav_lev',  len(z_interp) )
ncvar=ncOUT.createVariable(varmod,'f',('nsub','nav_lev'))
ncvar[:]=(STD)
ncOUT.close()

