#! /usr/bin/python
# source /gpfs/work/IscrC_MYMEDBIO/COPERNICUS/sequence.sh

import sys

from maskload import *

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle
import pylab as pl

Tr=np.loadtxt('transect.dat', dtype=trans_dt,skiprows=1)

jpk, jpj,jpi = tmask.shape
#LIST=[];
#StoredPointsForPlot=zeros(0,6);)

Matrices=[]
Minv=[]
LIST=np.array([],dtype=np.int)
FLUX_MAP=np.zeros((jpj,jpi), np.float32)
FLUX_MAP[~tmask[0,:,:]]=np.nan
nTr = len(Tr) if Tr.shape else 1
for iTrans in range(nTr):
    T = Tr[iTrans].copy() if Tr.shape else Tr.copy()

    DIST=(glamt - T['lon1'])**2 + (gphit - T['lat1'])**2; indStart = np.array(np.unravel_index(np.argmin(DIST),DIST.shape),ind2d)
    DIST=(glamt - T['lon2'])**2 + (gphit - T['lat2'])**2; ind__End = np.array(np.unravel_index(np.argmin(DIST),DIST.shape),ind2d)
    
    
#    %%    %%%%%%%%%%%%%%%% longitudinal %%%%%%%%%%%%%%%%%%%%%%%
    if T['type'] == 'Meridional':
        
        if np.isnan(T['depth1']): T['depth1'] = 0.
        
        if np.isnan(T['depth2']): T['depth2'] = 6000.  

        indDepStart = getDepthIndex(nav_lev, T['depth1']);
        indDep__End = getDepthIndex(nav_lev, T['depth2']);
        
        nY = ind__End['lat'] - indStart['lat'] +1
        nZ = indDep__End - indDepStart +1
        indT = np.zeros((nZ,nY),int);

        print Tr[iTrans]['name']
        print( indStart['lon']) 
        print( ind__End['lon']) 
        print( indStart['lat']) 
        print( ind__End['lat']) 
        print( '-------------') 

#        % ---------------------------------------------------
        if indStart['lon'] == ind__End['lon']:
            pos = indStart['lon'];
        else:
            candidates=np.arange(indStart['lon'],ind__End['lon']+1)
            nc=ind__End['lon']-indStart['lon']+1
            nPoints=np.zeros(nc)
            
            for iLon in range(nc):
                a = tmask[:,:,indStart['lon']:ind__End['lon']+1].copy()
                nPoints[iLon]=a.sum()
            
            isorted =nPoints.argsort()
            pos=candidates[isorted[nc-2]];
#        % -------------------------------------------------

        for j in range(indStart['lat'],ind__End['lat']+1):

            FLUX_MAP[j,pos] = FLUX_MAP[j,pos] + iTrans+1
            for k in range(indDepStart,indDep__End+1):
                
                indT[k-indDepStart, j-indStart['lat']] = index[0,k,j,pos]
        
    
#    %%    %%%%%%%%%%%%%%%%%%%%%%%%%%% latitudinal %%%%%%%%%%%%%%%
    if T['type'] == 'Zonal':
        if np.isnan(T['depth1']): T['depth1'] = 0.
        
        if np.isnan(T['depth2']): T['depth2'] = 6000.

        indDepStart = getDepthIndex(nav_lev, T['depth1']);
        indDep__End = getDepthIndex(nav_lev, T['depth2']);
        
        nX = ind__End['lon'] - indStart['lon'] +1
        nZ = indDep__End - indDepStart +1
        indT = np.zeros((nZ,nX),int);

        print Tr[iTrans]['name']
        print( indStart['lon'])
        print( ind__End['lon'])
        print( indStart['lat'])
        print( ind__End['lat'])
        print( '-------------')

        for i in range(indStart['lon'],ind__End['lon']+1):
            pos = indStart['lat'] #np.argmin(abs(gphiv[:,i] - T['lat1']));
            FLUX_MAP[pos,i]=FLUX_MAP[pos,i] + iTrans+1

            for k in range(indDepStart,indDep__End+1):
                indT[k-indDepStart, i-indStart['lon']] = index[0,k,pos,i];

#    %% horizontal plane
    if T['type'] == 'Vertical':
        nY=jpj
        nX=jpi
        k=getDepthIndex(nav_lev, T['depth1'])
        indT = np.zeros((nY,nX),int);
        indT = index[0,k,:,:]


#    if (T.depth1 == T.depth2)
#        T.type = 'horizontal';
#        T.Xriq = [T.lon1 T.lon2 T.lon2 T.lon1 T.lon1];
#        T.Yriq = [T.lat1 T.lat1 T.lat2 T.lat2 T.lat1];
#
#        ii = inpolygon(M.glamt,M.gphit,T.Xriq,T.Yriq);
#
#        indDepStart = getDepthIndex(M.nav_lev, T.depth1);
#        indT=squeeze(B.index(indDepStart,:,:));
#        indT(~ii) = 0;
#
#    end
#    %%

#    run Transect_plot_pl
#    disp('Type any key to continue'); pause
    ii=indT>0

    Matrices.append(indT.copy())
    indTr = indT[ii]
    LIST=np.append(LIST,indTr)


#FLUX_MAP[~tmask[0,:,:]]=np.nan
fig,ax=pl.subplots()
im=ax.pcolormesh(glamt,gphit, np.ma.masked_invalid(FLUX_MAP))
for iTrans in range(nTr):
    x=[Tr[iTrans]['lon1'],Tr[iTrans]['lon2']]
    y=[Tr[iTrans]['lat1'],Tr[iTrans]['lat2']]
    ax.plot(x,y,'r')
pl.colorbar(im)
fig.show()
LIST = np.unique(LIST)
#

for tr in range(nTr):
    indT = Matrices[tr]
    rINV = np.zeros(indT.shape)
    ii=   indT>0
    
    i,j=ii.nonzero()
    for t in range(len(i)):
        rINV[i[t],j[t]]=(LIST==indT[i[t],j[t]]).nonzero()[0]
    
    Minv.append(rINV.copy())

# writing the files

fid = open('Matrices.pkl','wb')
pickle.dump(Matrices,fid)
fid.close()

fid = open('rINV.pkl','wb')
pickle.dump(rINV,fid)
fid.close()

n=len(LIST)

print n

ncOUT=NC.netcdf_file('Fluxes.nc',"w")

ncOUT.createDimension('n',n)

ncvar=ncOUT.createVariable('index','i',('n',))

ncvar[:]=LIST

ncOUT.close()
