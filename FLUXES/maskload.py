
import V3C_mask as mask
import netCDF4 as NC4
import numpy as np

def getDepthIndex(nav_lev, lev_mask):
    jk_m = 0
    for jk in range(nav_lev.__len__()):
        if nav_lev[jk] < lev_mask:
            jk_m = jk
    return jk_m

#M=ncread(KB_mask.maskfile,'glamu','glamt','gphiv','gphit','nav_lev','tmask','umask','vmask');



M=NC4.Dataset(mask.maskfile,"r")
jpi=M.dimensions['x']
jpj=M.dimensions['y']
jpk=M.dimensions['z']


nav_lev =  M.variables['nav_lev'][:]
glamu     =  M.variables['glamu'][0,0,:,:]
glamt     =  M.variables['glamt'][0,0,:,:]
gphiv     =  M.variables['gphiv'][0,0,:,:]
gphit     =  M.variables['gphit'][0,0,:,:]
tmask   = (M.variables['tmask'][0,:,:,:]).astype(np.bool)
umask   = (M.variables['umask'][0,:,:,:]).astype(np.bool)
vmask   = (M.variables['vmask'][0,:,:,:]).astype(np.bool)


M.close()


B=NC4.Dataset(mask.bounmask,"r")

index =  B.variables['index'][:,:,:,:]
index_inv =  B.variables['index_inv'][:,:]
B.close()

print('done')
