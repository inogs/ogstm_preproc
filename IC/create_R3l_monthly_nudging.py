import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
#bitsea stuff
from bitsea.basins import V2 as OGS
from bitsea.commons.grid import RegularGrid
from bitsea.commons.mask import Mask
from IC import RSTwriter
from bitsea.commons import netcdf4
from bitsea.commons.interpolators import SeaOverLand

'''
create monthly reconstructed R3l fields from 
monthly satellite RRS412 files
according to a very shady algorythm
'''

def get_submask(mysub, lat, lon, tmask):
    seamask = tmask
    #seamask = ~np.isnan(X)
    tmask = mysub.is_inside(lon[:,np.newaxis], lat).T
    nanmask = tmask.astype(float)
    nanmask[~tmask] = np.nan
    area = RegularGrid(lon=lon, lat=lat).area
    # mask_land
    tmask = tmask * seamask
    nanmask = nanmask * seamask
    area = area * seamask * nanmask
    return tmask, nanmask, area

def saturate(Y, low_val=1.0, hig_val=1.0):
    Y = np.maximum(Y, low_val)
    Y = np.minimum(Y, hig_val)
    return Y

def CDOM_profile(surface_value, bottom_value, depth, depth0):
    k = 0.05
    res = surface_value + (bottom_value - surface_value) / (1.0 + np.exp(-k * (depth - depth0)))
    return res

#infile_nudge = '/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022_v17/R3l_yyyy0630-00:00:00.nc'
infile_nudge = '/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022_v20/R3l_yyyy0630-00:00:00.nc'

INDIR = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/SAT/'
infile_rrs_mean = 'rrs412_mediterranean_cmems_1999-2022.nc'
infile_rrs_mont_all = '????????_rrs_mediterranean_cmems_montly.nc'     #all moonths
infile_rrs_mont = '??0[789]????_rrs_mediterranean_cmems_montly.nc'  #summer (JAS) only
var_rrs = 'RRS412'

#OUTDIR = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/R3l_monthly/'
#OUTDIR = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/R3l_monthly-2107a/'
#OUTDIR = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/R3l_monthly-2107b/'
OUTDIR = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/R3l_monthly-2107e/'
#rstpfx='RST.19950101-00:00:00'
ndgpfx='R3l_yyyymm01-00:00:00.nc'

MASKFILE = '/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022_v6/wrkdir/MODEL/meshmask.nc'
TheMask = Mask.from_file(MASKFILE)
nav_lev = TheMask.zlevels
tmask = TheMask.mask.as_array()
tmask0 = tmask[0,:,:]
nanmask0 = tmask0.copy().astype(float)
nanmask0[~tmask0] = np.nan

#flist_rrs = sorted(glob(INDIR+infile_rrs_mont_all))
flist_rrs = sorted(glob(INDIR+infile_rrs_mont))
flist_rrs_all = sorted(glob(INDIR+infile_rrs_mont_all))

# compute mean RRS412
print('generating mean RRS412...')
D = xr.open_mfdataset(flist_rrs, combine='nested', concat_dim='time')
RRS_mean = D[var_rrs].mean(dim='time').values[:] #annual mean
D.close()

# loads R3l nudge file, takes the surface value
var_ndg = 'R3l'
D = xr.open_dataset(infile_nudge)
Lat_ndg = D['latitude'].values[:]
Lon_ndg = D['longitude'].values[:]
R3l = D[var_ndg].values[0,:,:]
D.close()

R3l_min = np.nanpercentile(R3l, 10) # =0.43
R3l_max = np.nanpercentile(R3l, 90) # =0.88
#R3l_5pc = np.nanpercentile(R3l, 5)
#R3l_5pc = np.nanmin(R3l)
#R3l_5pc = 0.15
R3l_5pc = 0.075
R3l_range = R3l_max - R3l_min

iRRS_mean= 1. / RRS_mean

tmask_alb, nanmask_alb, area_alb = get_submask(OGS.alb, Lat_ndg, Lon_ndg, tmask0)
tmask_atl, nanmask_atl, area_atl = get_submask(OGS.atl, Lat_ndg, Lon_ndg, tmask0)
tmask_nwm, nanmask_nwm, area_nwm = get_submask(OGS.nwm, Lat_ndg, Lon_ndg, tmask0)
tmask_lev3, nanmask_lev3, area_lev3 = get_submask(OGS.lev3, Lat_ndg, Lon_ndg, tmask0)

# compute regional means to decide iRRS range
very_hig_val = np.nansum(nanmask_alb * area_alb * iRRS_mean) / np.nansum(nanmask_alb * area_alb)
hig_val = np.nansum(nanmask_nwm * area_nwm * iRRS_mean) / np.nansum(nanmask_nwm * area_nwm)
low_val = np.nansum(nanmask_lev3 * area_lev3 * iRRS_mean) / np.nansum(nanmask_lev3 * area_lev3)

#iRRS = saturate(iRRS, low_val, hig_val)

iRRS_range = hig_val - low_val

# Compute reconstructed monthly R3l fields
# RESCALE iRRS so that it falls into the range of R3l
for ff in flist_rrs_all:
    print(ff)
    curr_year = ff.split('/')[-1][4:8]
    curr_month = ff.split('/')[-1][2:4]
    #outfile_rst = 
    outfile_ndg = ndgpfx.replace('yyyy',curr_year).replace('mm',curr_month)
    outfile_ndg = OUTDIR+outfile_ndg
    D = xr.open_dataset(ff)
    RRS = D[var_rrs].values[:]
    #
    iRRS = 1. / RRS
    # manually fill in Atlantic
    iRRS[tmask_atl[np.newaxis,:,:]] = very_hig_val
    iRRS = saturate(iRRS, low_val, very_hig_val)
    #X = ((iRRS - low_val) * (R3l_range / iRRS_range)) + np.nanmin(R3l)
    X = ((iRRS - low_val) * (R3l_range / iRRS_range)) + R3l_5pc
    X = SeaOverLand(X[0,:,:], 10)
    X = X * nanmask0
    # land points that should be water: 
    # II = np.isnan(X) * tmask0
    X_3D = CDOM_profile(X[np.newaxis,:,:], 1.0, nav_lev[:,np.newaxis,np.newaxis], 150.0)
    X_3D[tmask==0.0] = 1e20
    #RSTwriter(outfile_rst, "R3l", X_3D, TheMask)
    print(f'saving: '+outfile_ndg)
    netcdf4.write_3d_file(X_3D, "R3l", outfile_ndg, TheMask, compression=True)
    D.close()
    #

