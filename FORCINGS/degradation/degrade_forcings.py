import numpy as np
import xarray as xr
from glob import glob
import os
from pathlib import Path
from argparse import ArgumentParser

from itertools import chain
# my stuff
import degrade_mesh as dm
from commons import dump_netcdf, load_parameters
from bitsea.commons.mask import Mask
from previous import forcingswriter as FW

'''
degrades resolution of OGSTM physics forcings
on T, U, V, W grids
needs both original and degraded resolution meshmask.nc
(run degrade_mesh.py before)
input parameters (like directories etc.) go into a .yaml file
usage: python degrade_forcings.py -y degrade_forcings.yaml

yaml file example:
maskfile: '/path/to/meshmask_in.nc' #original mesh
maskfile_d: '/path/to/coarse/meshmask_out.nc' #degraded mesh
ffdir: '/path/to/wrkdir/MODEL/FORCINGS/'
outdir: '/path/to/output/'
ndeg: 6   #nr of cells to join in i, j
y0: 2000  #start year
yE: 2020  #end year

(takes about 3min per year of daily files on 1 node (256GB) and 16 cores,
if it doesn't die by OOM, which sometimes does and sometimes not. 
If it does go with 8 cores instead)
'''

def argument():
    parser = ArgumentParser()
    parser.add_argument('--yamlfile','-y',
                        type=str,
                        required=True,
                        help='file with all parameters')
    return parser.parse_args()


def load_mesh_light(maskfile, ndeg=1):
    '''
    load meshmask, expand fields
    just those needed to degrade forcings
    '''
    M = xr.open_dataset(maskfile)
    M1 = {}
    #
    M1["e1t"] = dm.xpnd_wrap(M["e1t"], 'interp', ndeg)
    M1["e1u"] = dm.xpnd_wrap(M["e1u"], 'interp', ndeg)
    M1["e1v"] = dm.xpnd_wrap(M["e1v"], 'interp', ndeg)
    M1["e2t"] = dm.xpnd_wrap(M["e2t"], 'edge', ndeg)
    M1["e2u"] = dm.xpnd_wrap(M["e2u"], 'edge', ndeg)
    M1["e2v"] = dm.xpnd_wrap(M["e2v"], 'edge', ndeg)
    M1["e3t_0"] = dm.xpnd_wrap(M["e3t_0"], 'edge', ndeg)
    M1["e3u_0"] = dm.xpnd_wrap(M["e3u_0"], 'edge', ndeg)
    M1["e3v_0"] = dm.xpnd_wrap(M["e3v_0"], 'edge', ndeg)
    M1["glamt"] = dm.xpnd_wrap(M["glamt"], 'interp', ndeg)
    M1["glamu"] = dm.xpnd_wrap(M["glamu"], 'interp', ndeg)
    M1["glamv"] = dm.xpnd_wrap(M["glamv"], 'interp', ndeg)
    M1["gphit"] = dm.xpnd_wrap(M["gphit"], 'interp', ndeg)
    M1["gphiu"] = dm.xpnd_wrap(M["gphiu"], 'interp', ndeg)
    M1["gphiv"] = dm.xpnd_wrap(M["gphiv"], 'interp', ndeg)
    M1["tmask"] = dm.xpnd_wrap(M["tmask"], 'edge', ndeg) #there is some error here!
    M1["umask"] = dm.xpnd_wrap(M["umask"], 'edge', ndeg) #there is some error here!
    M1["vmask"] = dm.xpnd_wrap(M["vmask"], 'edge', ndeg) #there is some error here!
    # Bonus variable for degrading forcings
    M1['h_column_t'] = M1['e3t_0'].sum(dim='z')
    #
    M.close()
    M1 = xr.Dataset(M1)
    return M1


def load_tfile(infile, ndeg=1):
    '''
    loads t-grid forcing file, pads a rim all around it
    to make j, i multiples of ndeg,
    populates rim with edge values
    (ugly, but it's the atlantic buffer that gets populated,
    the model won't use those values anyway)
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    # most of these ogstm does not need (commented)
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['votemper'] = dm.xpnd_wrap(F['votemper'], 'edge', ndeg)
    F1['vosaline'] = dm.xpnd_wrap(F['vosaline'], 'edge', ndeg)
    F1['sossheig'] = dm.xpnd_wrap(F['sossheig'], 'edge', ndeg)
    #F1['sossh_ib'] = dm.xpnd_wrap(F['sossh_ib'], 'edge', ndeg)
    #F1['sowaflup'] = dm.xpnd_wrap(F['sowaflup'], 'edge', ndeg)
    #F1['soevapor'] = dm.xpnd_wrap(F['soevapor'], 'edge', ndeg)
    #F1['soprecip'] = dm.xpnd_wrap(F['soprecip'], 'edge', ndeg)
    F1['sorunoff'] = dm.xpnd_wrap(F['sorunoff'], 'edge', ndeg)
    F1['soshfldo'] = dm.xpnd_wrap(F['soshfldo'], 'edge', ndeg)
    #F1['sohefldo'] = dm.xpnd_wrap(F['sohefldo'], 'edge', ndeg)
    #F1['solofldo'] = dm.xpnd_wrap(F['solofldo'], 'edge', ndeg)
    #F1['sosefldo'] = dm.xpnd_wrap(F['sosefldo'], 'edge', ndeg)
    #F1['solafldo'] = dm.xpnd_wrap(F['solafldo'], 'edge', ndeg)
    F1['somxl010'] = dm.xpnd_wrap(F['somxl010'], 'edge', ndeg)
    # overwrite coordinates, else it complains
    for vv in F1.keys():
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    #
    F1['deptht'] = F['deptht']
    F1['deptht_bounds'] = F['deptht_bounds']
    F1['time_instant'] = F['time_instant']
    F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_ufile(infile, ndeg=1):
    '''
    as in load_tfile()
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vozocrtx'] = dm.xpnd_wrap(F['vozocrtx'], 'edge', ndeg)
    F1['sozotaux'] = dm.xpnd_wrap(F['sozotaux'], 'edge', ndeg)
    vlist = list(F1.keys())
    # overwrite coordinates, else it complains
    for vv in vlist:
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthu'] = F['depthu'] 
    F1['depthu_bounds'] = F['depthu_bounds'] 
    F1['time_instant'] = F['time_instant'] 
    F1['time_instant_bounds'] = F['time_instant_bounds'] 
    F1['time_counter'] = F['time_counter'] 
    F1['time_counter_bounds'] = F['time_counter_bounds'] 
    F1['time_centered'] = F['time_centered'] 
    F1['time_centered_bounds'] = F['time_centered_bounds'] 
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_vfile(infile, ndeg=1):
    '''
    as in load_tfile(), load_ufile
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vomecrty'] = dm.xpnd_wrap(F['vomecrty'], 'edge', ndeg)
    F1['sometauy'] = dm.xpnd_wrap(F['sometauy'], 'edge', ndeg)
    vlist = list(F1.keys())
    # overwrite coordinates, else it complains
    for vv in vlist:
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthv'] = F['depthv']
    F1['depthv_bounds'] = F['depthv_bounds']
    F1['time_instant'] = F['time_instant']
    F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1

def load_wfile(infile, ndeg=1):
    '''
    as in load_tfile(), load_ufile(), load vfile()
    '''
    F = xr.open_dataset(infile)
    F1 = {}
    #
    F1['nav_lat'] = dm.xpnd_wrap(F['nav_lat'], 'interp', ndeg)
    F1['nav_lon'] = dm.xpnd_wrap(F['nav_lon'], 'interp', ndeg)
    F1['vovecrtz'] = dm.xpnd_wrap(F['vovecrtz'], 'edge', ndeg)
    F1['votkeavt'] = dm.xpnd_wrap(F['votkeavt'], 'edge', ndeg)
    # overwrite coordinates, else it complains
    for vv in F1.keys():
        F1[vv].coords['nav_lat'].values[:] = F1['nav_lat'].values[:]
        F1[vv].coords['nav_lon'].values[:] = F1['nav_lon'].values[:]
    F1['depthw'] = F['depthw']
    F1['depthw_bounds'] = F['depthw_bounds']
    #F1['time_instant'] = F['time_instant']
    #F1['time_instant_bounds'] = F['time_instant_bounds']
    F1['time_counter'] = F['time_counter']
    F1['time_counter_bounds'] = F['time_counter_bounds']
    F1['time_centered'] = F['time_centered']
    F1['time_centered_bounds'] = F['time_centered_bounds']
    #
    F1 = xr.Dataset(F1)
    F.close()
    return F1


def get_weights(M, T):
    '''
    given the meshmask and the current t-grid file (with free surface)
    computes cell surface areas on U, V grids and cell volume on T grid
    and surface area on T / W grid
    algorithm from ogstm/src/IO/forcing_phys.f90
    (NB, this one pure numpy because xarray is not good at broadcasting)

    Arguments:
    M : xarray dataset output of load_mesh_light()
    T : xarray dataset output of load_tfile()

    '''
    jpk = M.dims['z']
    h_column_t = M['h_column_t'].values[:]
    tmask =           M['tmask'].values[:]
    umask =           M['umask'].values[:]
    vmask =           M['vmask'].values[:]
    e3t_0 =           M['e3t_0'].values[:]
    e3u_0 =           M['e3u_0'].values[:]
    e3v_0 =           M['e3v_0'].values[:]
    e1u =               M['e1u'].values[:]
    e2u =               M['e2u'].values[:]
    e1v =               M['e1v'].values[:]
    e2v =               M['e2v'].values[:]
    e1t =               M['e1t'].values[:]
    e2t =               M['e2t'].values[:]
    At = e1t * e2t #e1v * e2u
    Aw = np.repeat(At, repeats=tmask.shape[1], axis=1)

    tmask_0=tmask[0,0,:,:].astype(bool)

    ssh = tmask[0,0,:,:] * T['sossheig'].values[:] # (1, y, x)
    ssh[0,~tmask_0] = 0.0 # avoiding nans in land processors
    correction_e3t = 1.0 + (ssh / h_column_t) # (1, y, x)
    e3t = e3t_0 * correction_e3t # broadcasting, resulting in (1, z, y, x)

    # corresponds to:
    # for jk in range(jpk):
    #     e3t[:,jk,:,:] = e3t_0[0,jk,:,:] * correction_e3t[0, : , :]


    e1u_x_e2u = (e1u * e2u).repeat(jpk, axis=1) #(1,jpk,jpj,jpi)
    e1v_x_e2v = (e1v * e2v).repeat(jpk, axis=1)
    e1t_x_e2t = (e1t * e2t).repeat(jpk, axis=1)
    diff_e3t = e3t - e3t_0
    
    assert np.isnan(diff_e3t).any() == False, "NaN values found in diff_e3t!"
    
    s0 = e1t_x_e2t * diff_e3t # ((1, jpk, jpj, jpi))
    s1 = np.zeros(s0.shape)
    s2 = np.zeros(s0.shape)
    s1[:,:,:,0:-1] = e1t_x_e2t[:,:,:,1:] * diff_e3t[:,:,:,1:]
    s2[:,:,0:-1,:] = e1t_x_e2t[:,:,1:,:] * diff_e3t[:,:,1:,:]

    e3u = e3u_0 + (0.5 * (umask / e1u_x_e2u) * (s0 + s1))
    e3v = e3v_0 + (0.5 * (vmask / e1v_x_e2v) * (s0 + s2))

        #  DO ji = 1,jpim1
        #  DO jj = 1,jpjm1
        #  DO jk = 1,jpk
        #      s0= e1t_x_e2t(jj,ji ) * diff_e3t(jk,jj,ji)
        #      s1= e1t_x_e2t(jj,ji+1) * diff_e3t(jk,jj,ji+1)
        #      s2= e1t_x_e2t(jj+1,ji) * diff_e3t(jk,jj+1,ji)
        #      e3udta(jk,jj,ji,2) = 0.5*(umask(jk,jj,ji)/(e1u_x_e2u(jj,ji)) * (s0 + s1))
        #      e3vdta(jk,jj,ji,2) = 0.5*(vmask(jk,jj,ji)/(e1v_x_e2v(jj,ji)) * (s0 + s2))
        #  ENDDO
        #  ENDDO
        #  ENDDO
    #
    B = {}
    B['e1v'] = xr.DataArray(e1v, name='e1v' , dims=('time','z_a','y','x'))
    B['e2u'] = xr.DataArray(e2u, name='e2u' , dims=('time','z_a','y','x'))
    B['At'] = xr.DataArray(At, name='At' , dims=('time','z_a','y','x'))
    B['Aw'] = xr.DataArray(Aw, name='Aw' , dims=('time','z','y','x'))
    B['Au'] = xr.DataArray((e2u * e3u), name='Au' , dims=('time','z','y','x'))
    B['Av'] = xr.DataArray((e1v * e3v), name='Av' , dims=('time','z','y','x'))
    B['V'] = xr.DataArray((e1t * e2t * e3t), name='V' , dims=('time','z','y','x')) #for degrading T-grid, don't weight land values
    B = xr.Dataset(B)
    return B



def vwmean2d(X, tmask_in, W, Mask_out):
    '''
    W (time: 1, z: 141, y: 66, y_b: 6, x: 219, x_b: 6)>
    mean, weighted by surface area or volume(t-grid)
    '''

    X = X.where(tmask_in.isel(z=0)) # set land points to nan
    weighted_sum = (W * X).sum(dim=('y_b', 'x_b'), skipna=True)
    weight_sum = W.sum(dim=('y_b', 'x_b'), skipna=True)
    weight_sum = weight_sum.where(Mask_out.mask[0,:], other=1.0)

    wmean = weighted_sum / weight_sum
    wmean = wmean.where(Mask_out.mask[0,:], other=1.e+20)
    assert np.isnan(wmean).any().values == False, "NaN values found in weighted mean!"

    return wmean

def vwmean(X, tmask_in, W, Mask_out):
    '''
    W (time: 1, z: 141, y: 66, y_b: 6, x: 219, x_b: 6)>
    mean, weighted by surface area or volume(t-grid)
    '''
    X = X.where(tmask_in) # set land points to nan
    weighted_sum = (W * X).sum(dim=('y_b', 'x_b'), skipna=True)
    weight_sum = W.sum(dim=('y_b', 'x_b'), skipna=True)
    weight_sum = weight_sum.where(Mask_out.mask, other=1.0) # to avoid division by zero
    
    wmean = weighted_sum / weight_sum
    wmean = wmean.where(Mask_out.mask, other=1.e+20)
    assert np.isnan(wmean).any().values == False, "NaN values found in weighted mean!"

    return wmean

    # OUT = xr.ones_like(X.isel(x_b=0, y_b=0)) * 1.e+20
    # _, jpk, jpj, jpi = OUT.shape
    # for ji in range(jpi):
    #     print(ji)
    #     for jj in range(jpj):
    #         for jk in range(jpk):
    #             if Mask_out.mask[jk,jj,ji]:
    #                 s = X.isel(z=jk, y=jj, x=ji).values
    #                 w = W.isel(z=jk, y=jj, x=ji).values
    #                 mask = tmask_in.isel(z=jk, y=jj, x=ji).values
    #                 s = s[mask]
    #                 w= w[mask]
    #                 OUT.values[0,jk,jj,ji] = (s * w).sum() / w.sum()

    #return OUT

def coarsen_U(U, umask_in, Au, Mask_out:Mask):
    '''
    Coarsens U field by weighted averaging the eastern column of each block of
    fine U-gridded array.

    Arguments:
    U        : xr.DataArray of fine U, reshaped in blocks (time, z, y, x, y_b, x_b)
    umask_in : xr.DataArray, reshaped in blocks umask of original mesh
    Au       : xr.DataArray, reshaped in blocks,
                with no nans,
                side U area, (e2u * e3u) in get_weights()

    Mask_out : Mask object, degraded mesh mask, having umask in 'mask' field

    Returns:
    U_coarse : xr.DataArray of coarse U, (time, z, y, x)
    '''

    U = U.where(umask_in) # set u-grid points to nan
    # take the eastern column of each block (U side)
    Au_y_b = Au.isel({'x_b':-1})
    U_y_b = U.isel({'x_b':-1})

    # Here skipna=True is important to exclude umask=0 points
    Flux_fine = (Au_y_b * U_y_b).sum(dim='y_b', skipna=True)

    # Here we take in account also the area of the umask=0 points
    Area_coarse_cells = Au_y_b.sum(dim='y_b')

    U_coarse = Flux_fine / Area_coarse_cells

    U_coarse = U_coarse.where(Mask_out.mask, other=1.e+20)
    assert np.isnan(U_coarse).any().values == False, "NaN values found in U coarse!"
    return U_coarse

def coarsen_V(V, vmask_in, Av, Mask_out:Mask):
    '''
    Coarsens V field by weighted averaging the northern column of each block of
    fine V-gridded array.

    Arguments:
    U        : xr.DataArray of fine U, reshaped in blocks (time, z, y, x, y_b, x_b)
    umask_in : xr.DataArray, reshaped in blocks umask of original mesh
    Au       : xr.DataArray, reshaped in blocks,
                with no nans,
                side U area, (e2u * e3u) in get_weights()

    Mask_out : Mask object, degraded mesh mask, having umask in 'mask' field

    Returns:
    U_coarse : xr.DataArray of coarse U, (time, z, y, x)
    '''

    V = V.where(vmask_in) # set v-grid points to nan
    # take the northern column of each block (V side)
    Av_y_b = Av.isel({'y_b':-1})
    V_y_b = V.isel({'y_b':-1})

    # Here skipna=True is important to exclude vmask=0 points
    Flux_fine = (Av_y_b * V_y_b).sum(dim='x_b', skipna=True)

    # Here we take in account also the area of the vmask=0 points
    Area_coarse_cells = Av_y_b.sum(dim='x_b')

    V_coarse = Flux_fine / Area_coarse_cells

    V_coarse = V_coarse.where(Mask_out.mask, other=1.e+20)
    assert np.isnan(V_coarse).any().values == False, "NaN values found in V coarse!"
    return V_coarse

def coarsen_taux(tau, umask_in, e2u, Mask_out):
    '''
    Coarsens taux field by weighted averaging the eastern column of each block of
    fine U-gridded array.

    Arguments:
    tau.     : xr.DataArray, reshaped in blocks (1, z, y, x, y_b, x_b)
    umask_in : xr.DataArray, reshaped in blocks umask of original mesh
    e2u      : xr.DataArray, reshaped in blocks
    Mask_out : Mask object, degraded mesh mask, having umask in 'mask' field
    Returns:
    tau_coarse : xr.DataArray of coarse tau, (1, 1, y, x)
    '''
    tau = tau.where(umask_in.isel(z=0)) # set u-grid points to nan
    # take the eastern column of each block (U side)
    e2u_y_b = e2u.isel({'x_b':-1})
    tau_y_b = tau.isel({'x_b':-1})

    # Here skipna=True is important to exclude umask=0 points
    Force_fine = (e2u_y_b * tau_y_b).sum(dim='y_b', skipna=True)

    # Here we take in account also the area of the umask=0 points
    Area_coarse_cells = e2u_y_b.sum(dim='y_b')
    tau_coarse = Force_fine / Area_coarse_cells

    tau_coarse = tau_coarse.where(Mask_out.mask[0,:], other=1.e+20)
    assert np.isnan(tau_coarse).any().values == False, "NaN values found in tau coarse!"
    return tau_coarse

def coarsen_tauy(tau, vmask_in, e1v, Mask_out):
    '''
    Coarsens tauy field by weighted averaging the eastern column of each block of 
    fine V-gridded array.

    Arguments:
    tau.      : xr.DataArray, reshaped in blocks (1, z, y, x, y_b, x_b)
    vmask_in : xr.DataArray, reshaped in blocks vmask of original mesh
    e1v      : xr.DataArray, reshaped in blocks
    Mask_out : Mask object, degraded mesh mask, having umask in 'mask' field
    Returns:
    tau_coarse : xr.DataArray of coarse tau, (1, 1, y, x)
    '''
    tau = tau.where(vmask_in.isel(z=0)) # set v-grid points to nan
    # take the eastern column of each block (V side)
    e1v_y_b = e1v.isel({'y_b':-1})
    tau_y_b = tau.isel({'y_b':-1})

    # Here skipna=True is important to exclude umask=0 points
    Force_fine = (e1v_y_b * tau_y_b).sum(dim='x_b', skipna=True)
    # Here we take in account also the area of the umask=0 points
    Area_coarse_cells = e1v_y_b.sum(dim='x_b')
    tau_coarse = Force_fine / Area_coarse_cells

    tau_coarse = tau_coarse.where(Mask_out.mask[0,:], other=1.e+20)
    assert np.isnan(tau_coarse).any().values == False, "NaN values found in tau coarse!"
    return tau_coarse


def degrade_V(V, vmask_in, Av, e1v, Mask_out, outfile, ndeg=1):
    '''
    degrades resolution of V-grid file and dumps it to netcdf
    Arguments:
    V        : xr.Dataset, v-grid forcing file, expanded but not yet reshaped in blocks
    vmask_in : xr.DataArray, reshaped in blocks vmask of original mesh
    Av       : xr.DataArray, reshaped in blocks cell surface area on v-grid
               (e1v * e3v) in get_weights(), with no nans
    e1v      : xr.DataArray, reshaped in blocks side V area, with no nans
    Mask_out : Mask object, degraded mesh mask
    outfile  : str, path to output file
    ndeg     : int, degradation factor

    '''
    V = V.rename({'time_counter':'time', 'depthv':'z'})
    Vd = {}
    Vr = dm.reshape_blocks(V['vomecrty'], ndeg)
    taur = dm.reshape_blocks(V['sometauy'], ndeg)

    Vd['vomecrty'] = coarsen_V(Vr, vmask_in, Av, Mask_out)
    Vd['sometauy'] = coarsen_tauy(taur, vmask_in, e1v, Mask_out)

    FW.writefileV(outfile,
                  Mask_out,
                  Vd['vomecrty'].values[:],
                  Vd['sometauy'].values[:]
                  )

def degrade_U(U, umask_in, Warea, e2u, Mask_out, outfile, ndeg=1):
    '''
    degrades resolution of U-grid file and dumps it to netcdf
    Arguments:
    U        : xr.Dataset, u-grid forcing file, expanded but not yet reshaped in blocks
    umask_in : xr.DataArray, reshaped in blocks umask of original mesh
    Warea    : xr.DataArray, reshaped in blocks cell surface area on u-grid
               (e2u * e3u) in get_weights(), with no nans
    e2u      : xr.DataArray, reshaped in blocks side U area, with no nans
    Mask_out : Mask object, degraded mesh mask
    outfile  : str, path to output file
    ndeg     : int, degradation factor
    '''
    U = U.rename({'time_counter':'time', 'depthu':'z'})
    Ud = {}
    
    Ur = dm.reshape_blocks(U['vozocrtx'], ndeg)
    taur = dm.reshape_blocks(U['sozotaux'], ndeg)

    Ud['vozocrtx'] = coarsen_U(Ur, umask_in, Warea, Mask_out)
    Ud['sozotaux'] = coarsen_taux(taur, umask_in, e2u, Mask_out)

    FW.writefileU( outfile,
                  Mask_out,
                  Ud['vozocrtx'].values[:],
                  Ud['sozotaux'].values[:]
                  )
    
    
def degrade_W(W, tmask_in, Wa, Mask_out, outfile, ndeg=1):
    '''
    degrades resolution of W-grid file and dumps it to netcdf

    Arguments:
    tmask_in : xarray DataArray, reshaped in blocks tmask of original mesh
    Wa       : xarray DataArray, reshaped in blocks cell surface area on t-grid
    Mask_out : Mask object, degraded mesh mask
    outfile  : str, path to output file
    ndeg     : int, degradation factor

    '''
    W = W.rename({'time_counter':'time', 'depthw':'z'})
    Wd = {}

    for var3d in ["vovecrtz", "votkeavt"]:
        X = dm.reshape_blocks(W[var3d], ndeg)
        Wd[var3d] = vwmean(X, tmask_in, Wa, Mask_out)

    FW.writefileW(outfile, 
                  Mask_out,
                  Wd['vovecrtz'].values,
                  Wd['votkeavt'].values)
    


def degrade_T(T,tmask_in, Wa, Wv, Mask_out, outfile, ndeg=1):
    '''
    degrades resolution of T-grid file and dumps it to netcdf
    Arguments:
    T        : xarray Dataset, t-grid forcing file, expanded but not yet reshaped in blocks
    tmask_in : xarray DataArray, reshaped in blocks tmask of original mesh
    Wa       : xarray DataArray, reshaped in blocks cell surface area on t-grid
    Wv       : xarray DataArray, reshaped in blocks cell volume on t-grid
    Mask_out : Mask object, degraded mesh mask
    outfile  : str, path to output file
    ndeg     : int, degradation factor
    '''
    T = T.rename({'time_counter':'time', 'deptht':'z'})
    Td = {}

    for var3d in ["votemper", "vosaline"]:
        X = dm.reshape_blocks(T[var3d], ndeg)
        Td[var3d] = vwmean(X, tmask_in, Wv, Mask_out)
    
    for var2d in ["sossheig", "soshfldo", "sorunoff", "somxl010"]:
        X = dm.reshape_blocks(T[var2d], ndeg)
        Td[var2d] = vwmean2d(X, tmask_in, Wa, Mask_out)

    FW.writefileT(outfile, 
                  Mask_out,
                  Td['votemper'].values,
                  Td['vosaline'].values,
                  Td['sossheig'].values,
                  Td['soshfldo'].values)
    

def make_outdir(outdir, outfile):
    #/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/qDEG_SETUP/FORCINGS/
    #T20020308-12:00:00.nc
    yyyy = outfile[1:5]
    mm = outfile[5:7]
    outdir = f'{outdir}/{yyyy}/{mm}/'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True) 
    else:
        pass
    return outdir

    
    


def get_flist(tuvw:str, Params:dict):
    '''
    
    '''
    ffdir = Params['ffdir']
    y0 = Params['y0']
    yE = Params['yE']
    flist = [glob(f'{ffdir}/{YYYY}/??/{tuvw}*.nc') for YYYY in range(y0, yE+1)]
    flist = sorted(list(chain.from_iterable(flist)))
    return flist

if __name__=='__main__':
    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nranks = comm.size
    except:
        comm = None
        rank = 0
        nranks = 1
    #
    Params = load_parameters(argument().yamlfile)
    #
    maskfile = Path(Params['maskfile'])
    maskfile_d = Path(Params['maskfile_d'])
    ffdir =  Path(Params['ffdir'])
    outdir = Path(Params['outdir'])
    y0 = Params['y0']
    yE = Params['yE']
    #
    ndeg = Params['ndeg']
    #
    flistt = get_flist('T', Params) 
    flistu = get_flist('U', Params)
    flistv = get_flist('V', Params)
    flistw = get_flist('W', Params)


    #
    # load mesh and expand lat, lon to make them multiples of ndeg
    
    print('loading mask')
    M = load_mesh_light(maskfile, ndeg) #NB, here Au, Av, V, are calculated with e3tuv_0
    Mask_out = Mask.from_file(maskfile_d)
    Mask_out_u = Mask.from_file(maskfile_d, mask_var_name='umask')
    Mask_out_v = Mask.from_file(maskfile_d, mask_var_name='vmask')

    tmask_in= dm.reshape_blocks(M['tmask'].astype(bool), ndeg)
    umask_in= dm.reshape_blocks(M['umask'].astype(bool), ndeg)
    vmask_in= dm.reshape_blocks(M['vmask'].astype(bool), ndeg)
    #
    # LOOP ON FILES
    print('loading T, U, V files')
    ziter = list(zip(flistt, flistu, flistv, flistw))
    for tfile, ufile, vfile, wfile in ziter[rank::nranks]:
        V = load_vfile(vfile, ndeg)
        U = load_ufile(ufile, ndeg)
        W = load_wfile(wfile, ndeg)
        T = load_tfile(tfile, ndeg)


        Weight = get_weights(M, T)

        
        Warea = dm.reshape_blocks(Weight['At'], ndeg)
        Warea_w = dm.reshape_blocks(Weight['Aw'], ndeg)
        Wvolume = dm.reshape_blocks(Weight['V'], ndeg)
        Warea_u = dm.reshape_blocks(Weight['Au'], ndeg)
        Warea_v = dm.reshape_blocks(Weight['Av'], ndeg)
        e2u = dm.reshape_blocks(Weight['e2u'], ndeg)
        e1v = dm.reshape_blocks(Weight['e1v'], ndeg)

        print('degrado')
        outft = tfile.split('/')[-1]
        outfu = ufile.split('/')[-1]
        outfv = vfile.split('/')[-1]
        outfw = wfile.split('/')[-1]
        outdir = make_outdir(outdir, outft)
        
        degrade_T(T, tmask_in, Warea, Wvolume, Mask_out, outdir + outft, ndeg)
        degrade_U(U, umask_in, Warea_u, e2u, Mask_out_u, outdir + outfu, ndeg)
        degrade_V(V, vmask_in, Warea_v, e1v, Mask_out_v, outdir + outfv, ndeg)
        degrade_W(W, tmask_in, Warea_w, Mask_out, outdir + outfw, ndeg)


