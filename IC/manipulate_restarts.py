import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Performs these requests
    in ion1, ion2, ion3, lev1, lev2, lev3, lev4: N1p=N1p*0.75 in 50-500 metri  e   N1p=N1p*0.95 in 500-1000 metri
    in nwm N1p=N1p*1.2 0-500 e N1p=N1p*1.1 500-1000
    in nwm N3n=1 0-50; N3n=3 50-100; N3n=5 100-150; N3n=6.5 150-200; N3n=7.0 200-500


    O3h=O3h*0.98 in swm1 swm2 nwm tyr1 tyr2 nello strato 0-50m
    O3c=O3c*0.98 in swm1 swm2 nwm tyr1 tyr2 nello strato 0-50m
    O3c=O3c*0.95 in ion1 nello strato 0-100m
    O3h=O3h*0.95 in ion1 nello strato 0-100m
    O3c=O3c*0.98 in ion2 nello strato 0-100m
    O3h=O3h*0.98 in ion2 nello strato 0-100m

    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/some/path/with/STAT_PROFILES')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output directory'
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )

    return parser.parse_args()

args = argument()

import numpy as np
from commons.layer import Layer
import pylab as pl
from basins import V2
from basins.basin import ComposedBasin
from commons import netcdf4
from commons.mask import Mask
from commons.submask import SubMask
from IC import RSTwriter
from commons.utils import addsep
TheMask=Mask(args.maskfile)

INPUTDIR=addsep(args.inputdir) #"/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/RESTARTS/"
OUTDIR  =addsep(args.outdir) #"/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/"
jpk, jpj, jpi=TheMask.shape
N1p_to_modify=ComposedBasin('N1p',[V2.lev1,V2.lev2,V2.lev3,V2.lev4,V2.ion1, V2.ion2, V2.ion3],'reduction of N1p')
z=TheMask.zlevels
z2=z[:-1]+z[1:]

def get_profile_O3h(layer):
    values=np.ones((jpk,),np.float32)
    ii=(z<layer.bottom) & (z>layer.top)
    values[ii]=values[ii]*0.98
    for _ in range(5):
        v2=np.interp(z2, z, values, left=0.98)
        values=np.interp(z,z2,v2,left=0.98)
    values[:5]=0.98
    return values

def get_profile_N1p(layer):
    values=np.ones((jpk,),np.float32)
    ii=(z<layer.bottom) & (z>layer.top)
    values[ii]=values[ii]*0.8
    for _ in range(5):
        v2=np.interp(z2, z, values)
        values=np.interp(z,z2,v2)   
    return values

profile_O3hO3c =get_profile_O3h(Layer(0,20))
profile_N1p    =get_profile_N1p(Layer(50,500))





pl.close('all')
fig,ax=pl.subplots()
ax.plot(profile_O3hO3c,z,'r.-')
ax.plot(profile_N1p,z,'b.-')
ax.invert_yaxis()
fig.show()


TheMask_0= TheMask.cut_at_level(0)
S=SubMask(V2.ion1,maskobject=TheMask_0)

for var in ["O3c","O3h"]:
    inputfile=INPUTDIR + "RST.20161201-00:00:00." + var + ".nc"
    outfile  =  OUTDIR + "RST.20161201-00:00:00." + var + ".nc"
    R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]
    for i in range(jpi):
        for j in range(jpj):
            if S.mask[0,j,i]:
                R[:,j,i] = R[:,j,i]*profile_O3hO3c
    R[~TheMask.mask]=1.e+20
    RSTwriter(outfile, var, R, TheMask)

S=SubMask(N1p_to_modify,maskobject=TheMask_0)
for var in ["N1p"]:
    inputfile=INPUTDIR + "RST.20161201-00:00:00." + var + ".nc"
    outfile  =  OUTDIR + "RST.20161201-00:00:00." + var + ".nc"
    R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]
    for i in range(jpi):
        for j in range(jpj):
            if S.mask[0,j,i]:
                R[:,j,i] = R[:,j,i]*profile_N1p
    R[~TheMask.mask]=1.e+20
    RSTwriter(outfile, var, R, TheMask)    
    


