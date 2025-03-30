from bitsea.commons.Timelist import TimeInterval,TimeList
from bitsea.commons import season
import numpy as np
Seas_obj=season.season()
Seas_obj.setseasons(["0101","0701"], ["winter","summer"])


maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V4/meshmask.nc"
INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/REA/UNZIPPED/"
TI = TimeInterval("1999","2015","%Y")

maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc"
maskfile="/gpfs/scratch/userexternal/gbolzon0/TRANSITION_24/meshmask.nc"
INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/TRANSITION_24/AVE/"
TI = TimeInterval("2016","2018","%Y")

VARLIST=['N1p', 'N3n', 'N5s','O3c','O3h','HgII','MMHg']
mydtype=[(var,np.float32) for var in VARLIST]
