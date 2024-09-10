import glob,os
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
from IC import RSTwriter
import numpy as np

INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/REA_24/TEST_13/wrkdir/MODEL/AVE_FREQ_2"
RST_DIR= "/gpfs/scratch/userexternal/gbolzon0/REA_24/TEST_13/wrkdir/MODEL/RESTARTS/"
rstlist = glob.glob(RST_DIR + "RST.1997*")

varlist=[]
for filename in rstlist:
    basename = os.path.basename(filename)
    var=basename[22:25]
    varlist.append(var)
varlist.remove('R7c')

TheMask = Mask("/gpfs/scratch/userexternal/gbolzon0/REA_24/TEST_13/wrkdir/MODEL/meshmask.nc")

for var in varlist:
    print var
    outfile = RST_DIR +  "RST.19980601-00:00:00." + var + ".nc"
    filename = INPUTDIR + "/ave.19980516-12:00:00." + var + ".nc"
    V = DataExtractor(TheMask,filename,var).values.astype(np.float64)
    
    RSTwriter(outfile, var,V,TheMask)
    
