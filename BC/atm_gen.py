import numpy as np
from basins import V2
from commons.mask import Mask
from commons.submask import SubMask

filename="atmDep_16sub_openAndCoast_winAndSum.txt"
maskfile="/Users/gbolzon/eclipse-workspace/preproc/IC/meshmask.nc"
TheMask=Mask(maskfile)
Mask0 = TheMask.cut_at_le

mydtype=[('sub','U20')]
VARS=['NOw','NCw','NOs','NCs','POw','POs','PCw','PCs']
mydtype.extend( [(s,np.float32) for s in VARS] )
A=np.loadtxt(filename,dtype=mydtype, skiprows=2)

A[A['sub']=='alb']['NCw']

for isub, sub in V2.Pred:
    