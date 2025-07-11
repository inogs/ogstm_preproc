import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
#bitsea stuff
from bitsea.basins import V2 as OGS
from bitsea.commons.grid import RegularGrid

'''
creates a csv file with mean surface R3l in each subbasin,
with gradient estimated from satellite reflectance and
absolute value from our tuned R3l values in nwm and lev3.
these values are then fed to ANOTHER SCRIPT that extends at depth and
create the 3D nudging file


SPiani notes

for mysub in OGS.P:
    print(mysub)

MM = mysub.is_inside(Lon_rrs[:,np.newaxis], Lat_rrs)

R = RegularGrid(lon=Lon_rrs, lat=Lat_rrs)
Area = R.area
'''

def get_submask(mysub, lat, lon, X):
    seamask = ~np.isnan(X)
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

PL_dir = '/g100_work/OGS23_PRACE_IT/plazzari/COPERNICUS/ANALYZE_REFLECTANCE/'
CA_dir = '/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022_v18/wrkdir/POSTPROC/output/AVE_FREQ_3/RRS_DAILY/2003/'
CS_dir = '/g100_work/OGS23_PRACE_IT/csoto/rrs_data/V11C/SAT/MONTLY/montly_20*/'

#infile_nudge = '/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022_v17/R3l_yyyy0630-00:00:00.nc'
infile_nudge = '/g100_work/OGS_devC/camadio/Neccton_hindcast1999_2022_v17/R3l_yyyy0630-00:00:00.nc'
infile_rrs = '*_rrs_mediterranean_cmems_montly.nc'

flist_rrs = sorted(glob(CS_dir+infile_rrs)) #this also gets in 2000 and 2001
#flist_rrs = [ff for ff in flist_rrs if ]

# loads mean RRS412 field
var_rrs = 'RRS412'
D = xr.open_mfdataset(flist_rrs, combine='nested', concat_dim='time')
Lat_rrs = D['lat'].values[:]
Lon_rrs = D['lon'].values[:]
RRS = D[var_rrs].mean(dim='time').values[:]
D.close()

#RRS_range = np.nanmax(RRS) - np.nanmin(RRS)
RRS_range = np.nanpercentile(RRS, 99) - np.nanpercentile(RRS, 1)

# loads R3l nudge file, takes the surface value
var_ndg = 'R3l'
D = xr.open_dataset(infile_nudge)
Lat_ndg = D['latitude'].values[:]
Lon_ndg = D['longitude'].values[:]
R3l = D[var_ndg].values[0,:,:]
D.close()

R3l_min = np.nanpercentile(R3l, 10) # =0.43
R3l_max = np.nanpercentile(R3l, 90) # =0.88
#R3l_max = 1.0
R3l_range = R3l_max - R3l_min

#X = ((RRS - np.nanmin(RRS)) * (R3l_range / RRS_range)) + np.nanmin(R3l)

iRRS = 1. / RRS

#hig_val = np.nanmean(iRRS[500:1100, 1100:]) # circa nwm
#hig_val = np.nanmean(iRRS[400:800,:500])   # circa alb 
#low_val = np.nanmean(iRRS[:400, 2500:])     # circa lev3/4

tmask_alb, nanmask_alb, area_alb = get_submask(OGS.alb, Lat_rrs, Lon_rrs, RRS)
tmask_nwm, nanmask_nwm, area_nwm = get_submask(OGS.nwm, Lat_rrs, Lon_rrs, RRS)
tmask_lev3, nanmask_lev3, area_lev3 = get_submask(OGS.lev3, Lat_rrs, Lon_rrs, RRS)

'''
either 
R3l_max = 1.0 & hig_val = hig_val(alb)
or
R3l_max = 0.88 % hig_val = hig_val(nwm)
'''

# compute regional means to decide iRRS range
#hig_val = np.nansum(nanmask_alb * area_alb * iRRS) / np.nansum(nanmask_alb * area_alb)
hig_val = np.nansum(nanmask_nwm * area_nwm * iRRS) / np.nansum(nanmask_nwm * area_nwm)
low_val = np.nansum(nanmask_lev3 * area_lev3 * iRRS) / np.nansum(nanmask_lev3 * area_lev3)

#iRRS = saturate(iRRS, low_val, hig_val)

iRRS_range = hig_val - low_val
#iRRS_range = np.nanmax(iRRS) - np.nanmin(iRRS)

# RESCALE iRRS so that it falls into the range of R3l
#X = ((iRRS - np.nanmin(iRRS)) * (R3l_range / iRRS_range)) + np.nanmin(R3l)
X = ((iRRS - low_val) * (R3l_range / iRRS_range)) + np.nanmin(R3l)

fig0 = plt.figure()
ax0 = plt.axes()
#im0 = ax0.pcolormesh(X, vmin=R3l_min, vmax=R3l_max)
im0 = ax0.pcolormesh(X, vmin=R3l_min, vmax=1.0)
plt.colorbar(im0)
fig0.suptitle('R3l new')
plt.savefig('R3l_from_iRRS.png')

fig1 = plt.figure()
ax1 = plt.axes()
im1 = ax1.pcolormesh(R3l, vmin=R3l_min, vmax=R3l_max)
plt.colorbar(im1)
fig1.suptitle('R3l old')
plt.savefig('R3l_2basin.png')

fig2 = plt.figure()
ax2 = plt.axes()
cmin = np.nanpercentile(iRRS, 10)
cmax = np.nanpercentile(iRRS, 90)
im2 = ax2.pcolormesh(iRRS, vmin=cmin, vmax=cmax)
plt.colorbar(im2)
fig2.suptitle('iRRS')
plt.savefig('iRRS.png')

# basin mean values
blist = [ii.name for ii in OGS.P.basin_list]
Vals_Dict = dict.fromkeys(blist)
for mysub in OGS.P:
    sn = mysub.name
    tmask, nanmask, area = get_submask(mysub, Lat_rrs, Lon_rrs, RRS)
    sub_mean = np.nansum(nanmask * area * X) / np.nansum(nanmask * area)
    print(f'{sn}: {sub_mean}')
    Vals_Dict[sn] = sub_mean

#manually put atl in
Vals_Dict['atl'] = Vals_Dict['alb']

#manually set aeg value because else it is too high!
Vals_Dict['aeg'] = 0.55

C = pd.DataFrame.from_dict(Vals_Dict, orient='index', columns=['R3l'])
outdir = '/g100_work/OGS23_PRACE_IT/ggalli00/NECCTON/SCRIPTS/'
outfile = outdir+'mean_R3l_from_RRS412.csv'
print('saving: '+outfile)
C.to_csv(outfile)


#for mysub in OGS.P:
#    sn = mysub.name
#    tmask, nanmask, area = get_submask(mysub, Lat_rrs, Lon_rrs, iRRS)
#    sub_mean = np.nansum(nanmask * area * iRRS) / np.nansum(nanmask * area)
#    print(f'{sn}: {sub_mean}')
