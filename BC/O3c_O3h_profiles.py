from pathlib import Path
from bitsea.timeseries.plot import read_pickle_file
from bitsea.commons.mask import Mask
from bitsea.commons.density import get_density_from_4D_profiles

maskfile=Path("/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask_CMCC.nc")
TheMask = Mask.from_file(maskfile, e3t_var_name="e3t_0")
fileS = Path("/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/vosaline.pkl")
fileT = Path("/g100_scratch/userexternal/gbolzon0/V13C/GIB_profiles/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/votemper.pkl")
TIMESERIES_T, TL = read_pickle_file(fileT)
TIMESERIES_S, _ = read_pickle_file(fileS)
rho = get_density_from_4D_profiles(TIMESERIES_T[:,:,:,:,0], TIMESERIES_S[:,:,:,:,0], TheMask)

iCoast = 2 # everywhere
rho=rho[:,:,iCoast,:]
rho_annual = rho.mean(axis=0)
