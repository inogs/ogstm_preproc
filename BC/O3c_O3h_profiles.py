from pathlib import Path
from bitsea.timeseries.plot import read_pickle_file
from bitsea.commons.mask import Mask
from bitsea.commons.density import get_density_from_4D_profiles
import numpy as np
import xarray as xr

maskfile=Path("/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc")
TheMask = Mask.from_file(maskfile, e3t_var_name="e3t_0")
maskfile_CMCC=Path("/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc")
CMCCMask = Mask.from_file(maskfile_CMCC, e3t_var_name="e3t_0")
jpk,_,_ = TheMask.shape
# T, S over the box Rectangle(-8.0,  -6.0,   35.0,   36.5)
fileS = Path("/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/vosaline.pkl")
fileT = Path("/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/votemper.pkl")
TIMESERIES_T, TL = read_pickle_file(fileT)
TIMESERIES_S, _ = read_pickle_file(fileS)
rho = get_density_from_4D_profiles(TIMESERIES_T[:,:,:,:,0], TIMESERIES_S[:,:,:,:,0], CMCCMask)

iCoast = 2 # everywhere
rho=rho[:,:,iCoast,:]
rho_annual = np.squeeze(rho.mean(axis=0))
rho_annual = np.asarray(rho_annual)
first_masked = np.nonzero(rho_annual==0)[0][0]
rho_annual[first_masked:] = rho_annual[first_masked-1]

# GLODAP Climatology profile lat35_36.2 lon-9_-6.5 # micromol Kg-1
glodap_depth = [12.5, 37.5, 75, 125, 175, 250, 350, 500, 700, 900, 1100, 1350, 1750, 2500]
glodap_tco2 = [2093.125, 2094.167, 2110.68, 2127.86, 2130.167, 2136.583, 2141.6, 
    2157.267, 2179.471, 2195.233, 2204.933, 2204.09, 2187.871, 2183.814]
glodap_talk = [2392.122, 2382.2, 2381.056, 2377.53, 2364.288, 2355.395, 2347.738, 
    2337.467, 2349.1, 2371.542, 2382.508, 2387.511, 2342.51, 2337.6 ] 

glodap_tco2_interp = np.interp(TheMask.zlevels, glodap_depth, glodap_tco2)
glodap_talk_interp = np.interp(TheMask.zlevels, glodap_depth, glodap_talk)

O3c = rho_annual[:jpk] * glodap_tco2_interp * 12 /1000
O3h = rho_annual[:jpk] * glodap_talk_interp /1000


dims = ["nav_lev"]
coords = {"nav_lev": TheMask.zlevels}


dataset = xr.Dataset(
    data_vars={
        "O3h": (tuple(dims), O3h.astype(np.float32), {"units": " [mmol m^-3]"}),
    },
    coords=coords,
)
dataset.to_netcdf(Path("O3h.nc"))
dataset = xr.Dataset(
    data_vars={
        "O3c": (tuple(dims), O3c.astype(np.float32), {"units": " [mmol m^-3]"}),
    },
    coords=coords,
)
dataset.to_netcdf(Path("O3c.nc"))