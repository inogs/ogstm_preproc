import numpy as np
from bitsea.basins import V2 as OGS
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.commons import netcdf4

filename="atmDep_16sub_openAndCoast_winAndSum.txt"
maskfile="/Users/gbolzon/eclipse-workspace/preproc/IC/meshmask.nc"
TheMask=Mask.from_file(maskfile)
jpk, jpj, jpi = TheMask.shape
Mask0 = TheMask.cut_at_level(0)
mask200= TheMask.mask_at_level(200)

mydtype=[('sub','U20')]
VARS=['NOw','NOs','NCw','NCs','POw','POs','PCw','PCs']
mydtype.extend( [(s,np.float32) for s in VARS] )
A=np.loadtxt(filename,dtype=mydtype, skiprows=2)

A[A['sub']=='alb']['NCw']

dtype = [(sub.name, bool) for sub in OGS.Pred]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.Pred:
    sbmask         = SubMask(sub, Mask0).mask
    SUB[sub.name]  = sbmask[0,:,:]

N3n_win = np.zeros((jpj,jpi),np.float32)
N3n_sum = np.zeros((jpj,jpi),np.float32)
N1p_win = np.zeros((jpj,jpi),np.float32)
N1p_sum = np.zeros((jpj,jpi),np.float32)


for sub in OGS.Pred:
    if sub.name=='atl':continue
    ii_sub = A['sub']==sub.name
    mask = SUB[sub.name]
    N3n_win[mask] = A[ii_sub]['NCw']
    N3n_sum[mask] = A[ii_sub]['NCs']
    N1p_win[mask] = A[ii_sub]['PCw']
    N1p_sum[mask] = A[ii_sub]['PCs']

    mask = SUB[sub.name] & mask200
    N3n_win[mask] = A[ii_sub]['NOw']
    N3n_sum[mask] = A[ii_sub]['NOs']
    N1p_win[mask] = A[ii_sub]['POw']
    N1p_sum[mask] = A[ii_sub]['POs']

netcdf4.write_2d_file(N3n_win, 'atm_N3n', 'ATM.yyyy0215-00:00:00.nc', TheMask, compression=True)
netcdf4.write_2d_file(N1p_win, 'atm_N1p', 'ATM.yyyy0215-00:00:00.nc', TheMask, compression=True)

netcdf4.write_2d_file(N3n_sum, 'atm_N3n', 'ATM.yyyy0815-00:00:00.nc', TheMask, compression=True)
netcdf4.write_2d_file(N1p_sum, 'atm_N1p', 'ATM.yyyy0815-00:00:00.nc', TheMask, compression=True)

