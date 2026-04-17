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
from bitsea.commons.layer import Layer
import pylab as pl
from bitsea.basins import V2
from bitsea.basins.basin import ComposedBasin
from bitsea.commons import netcdf4
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from IC import RSTwriter, smoother
from bitsea.commons.utils import addsep
TheMask=Mask.from_file(args.maskfile)

INPUTDIR=addsep(args.inputdir) #"/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/RESTARTS/"
OUTDIR  =addsep(args.outdir) #"/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_03/wrkdir/MODEL/"
jpk, jpj, jpi=TheMask.shape
z=TheMask.zlevels
z2=z[:-1]+z[1:]
def smooth(values):
    for _ in range(2):
        v2=np.interp(z2, z, values)
        values=np.interp(z,z2,v2)
    return values

def get_profile_O3h(layer, value):
    values=np.ones((jpk,),np.float32)
    ii=(z<layer.bottom) & (z>layer.top)
    values[ii]=values[ii]*value
    for _ in range(2):
        v2=np.interp(z2, z, values, left=value)
        values=np.interp(z,z2,v2,left=value)
        values[:3]=value
    return values


def get_profile1_N1p():
    L1=Layer(50,500)
    L2=Layer(500,1000)
    values=np.ones((jpk,),np.float32)
    ii=(z<L1.bottom) & (z>L1.top); values[ii]=values[ii]*0.75
    ii=(z<L2.bottom) & (z>L2.top); values[ii]=values[ii]*0.95
    return smooth(values)

def get_profile2_N1p():
    L1=Layer(0,500)
    L2=Layer(500,1000)
    values=np.ones((jpk,),np.float32)
    ii=(z<L1.bottom) & (z>L1.top); values[ii]=values[ii]*1.2
    ii=(z<L2.bottom) & (z>L2.top); values[ii]=values[ii]*1.1
    return smooth(values)

def get_profile_N3n():
    '''
        in nwm
        N3n=1 0-50
        N3n=3 50-100
        N3n=5 100-150
        N3n=6.5 150-200
        N3n=7.0 200-500
    '''
    L1=Layer(0,50)
    L2=Layer(50,100)
    L3=Layer(100,150)
    L4=Layer(150,200)
    L5=Layer(200,5000)
    values=np.ones((jpk,),np.float32)
    ii=(z<L1.bottom) & (z>L1.top); values[ii]=1
    ii=(z<L2.bottom) & (z>L2.top); values[ii]=3
    ii=(z<L3.bottom) & (z>L3.top); values[ii]=5
    ii=(z<L4.bottom) & (z>L4.top); values[ii]=6.5
    ii=(z<L5.bottom) & (z>L5.top); values[ii]=7.0
    for _ in range(1):
        v2=np.interp(z2, z, values)
        values=np.interp(z,z2,v2)   
    return values


profile1_O3hO3c =get_profile_O3h(Layer(0, 50),0.98) # sub O3h_to_modify
profile2_O3hO3c =get_profile_O3h(Layer(0,100),0.95) # on ion1
profile3_O3hO3c =get_profile_O3h(Layer(0,100),0.98) # on ion2

profile1_N1p    =get_profile1_N1p() #N1p_to_modify
profile2_N1p    =get_profile2_N1p() # on nwm
profile_N3n     = get_profile_N3n() # on nwm



pl.close('all')
fig,ax=pl.subplots()
ax.plot(profile1_O3hO3c,z,'r.-')
ax.plot(profile2_O3hO3c,z,'g.-')
ax.plot(profile3_O3hO3c,z,'b.-')

#ax.plot(profile1_N1p,z,'b.-')
#ax.plot(profile2_N1p,z,'r.-')
#ax.plot(profile_N3n,z,'b.-')
ax.set_ylim([0,500])
ax.invert_yaxis()
fig.show()



TheMask_0= TheMask.cut_at_level(0)
N1p_to_modify=ComposedBasin('N1p',[V2.lev1,V2.lev2,V2.lev3,V2.lev4,V2.ion1, V2.ion2, V2.ion3],'reduction of N1p')
O3h_to_modify=ComposedBasin('O3h',[V2.swm1,V2.swm2, V2.nwm, V2.tyr1 ],'reduction of O3h, O3h in 0-50')

S_O3h_1=SubMask(O3h_to_modify, TheMask_0)
S_O3h_2=SubMask(V2.ion1, TheMask_0)
S_O3h_3=SubMask(V2.ion2, TheMask_0)

datestr="20161001-00:00:00"

for var in ["O3c","O3h"]:
    inputfile="%sRST.%s.%s.nc" %(INPUTDIR,datestr,var)
    outfile  ="%sRST.%s.%s.nc" %(OUTDIR  ,datestr,var)
    R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]

    for i in range(jpi):
        for j in range(jpj):
            if S_O3h_1.mask[0,j,i]:
                R[:,j,i] = R[:,j,i]*profile1_O3hO3c
    for i in range(jpi):
        for j in range(jpj):
            if S_O3h_2.mask[0,j,i]:
                R[:,j,i] = R[:,j,i]*profile2_O3hO3c
    for i in range(jpi):
        for j in range(jpj):
            if S_O3h_3.mask[0,j,i]:
                R[:,j,i] = R[:,j,i]*profile3_O3hO3c
    RST_s = smoother(TheMask, R)
    RSTwriter(outfile, var, RST_s, TheMask)

var="N1p"
S_N1p_1=SubMask(N1p_to_modify, TheMask_0)
S_N1p_2=SubMask(V2.nwm, TheMask_0)

inputfile="%sRST.%s.%s.nc" %(INPUTDIR,datestr,var)
outfile  ="%sRST.%s.%s.nc" %(OUTDIR  ,datestr,var)
R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]
for i in range(jpi):
    for j in range(jpj):
        if S_N1p_1.mask[0,j,i]:
            R[:,j,i] = R[:,j,i]*profile1_N1p
for i in range(jpi):
    for j in range(jpj):
        if S_N1p_2.mask[0,j,i]:
            R[:,j,i] = R[:,j,i]*profile2_N1p
RST_s = smoother(TheMask, R)
RSTwriter(outfile, var, RST_s, TheMask)

var="N3n"
S_N3n=SubMask(V2.nwm, TheMask_0)

inputfile="%sRST.%s.%s.nc" %(INPUTDIR,datestr,var)
outfile  ="%sRST.%s.%s.nc" %(OUTDIR  ,datestr,var)
R=netcdf4.readfile(inputfile, "TRN" + var)[0,:,:,:]
for i in range(jpi):
    for j in range(jpj):
        if S_N3n.mask[0,j,i]:
            R[:,j,i] = profile_N3n
RST_s = smoother(TheMask, R)
RSTwriter(outfile, var, RST_s, TheMask)
    


