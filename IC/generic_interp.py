from bitsea.commons import netcdf4
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons import interpolators
import numpy as np
from IC import RSTwriter


Mask16 = Mask.from_file("/g100_scratch/userexternal/grosati0/SIM_ATM/SIM_newATM/wrkdir/MODEL/meshmask.nc")
Mask24 = Mask.from_file("/g100_work/OGS23_PRACE_IT/grosati/NECCTON/PREPROC_ogstm/ogstm_preproc/BC/meshmask.nc")


filename = (
"/g100_scratch/userexternal/grosati0/SIM_ATM/SIM_newATM/wrkdir/MODEL/AVE_FREQ_1/ave.20150620-00:00:00.N1p.nc"
)
VAR = DataExtractor(Mask16, filename, "N1p").values

A = interpolators.space_interpolator_griddata(Mask24, Mask16, VAR).astype(np.float64)
RSTwriter("RST.20220101-00:00:00.N1p.nc", "N1p", A, Mask24)
