# This script prints as standart output two columns of indexes
# that have to be putted in rivers xls file for reduced mesh

import numpy as np
class indexer():
    def __init__(self,I,J):
        self.I = I
        self.J = J
        self.jpi = len(I)
        self.jpj = len(J)
    def switch_lon_index(self,ifine):
        for k in range(self.jpi):
            if self.I[k]>ifine:
                coarse_ind = k
                break
        return coarse_ind-1
    def switch_lat_index(self,jfine):
        for k in range(self.jpj):
            if self.J[k]>jfine:
                coarse_ind = k
                break
        return coarse_ind-1    
     




I_ind = np.loadtxt('ogstm_South_West_I_indexes.txt',np.int32)
J_ind = np.loadtxt('ogstm_South_West_J_indexes.txt',np.int32)
Indexer = indexer(I_ind,J_ind)


# Prese le colonne I,J del file river_8_on_meshgrid_v2.xlsx e messe in nei files I_river_16_from_xls.txt  J_river_16_from_xls.txt
I_river_16 = np.loadtxt("I_river_16_from_xls.txt",np.int32) -1 + 222 # because plazzari removed 222
J_river_16 = np.loadtxt("J_river_16_from_xls.txt",np.int32) -1


riverPoints = len(I_river_16)
I_red = np.zeros((riverPoints,),np.int32)
J_red = np.zeros((riverPoints,),np.int32)

for iR in range(riverPoints):
    i = Indexer.switch_lon_index(I_river_16[iR])
    j = Indexer.switch_lat_index(J_river_16[iR])
    print "%d\t%d" %(i+1,j+1)
    I_red[iR] = i
    J_red[iR] = j


from commons.mask import Mask
TheMask = Mask("meshmask_469.nc")

from layer_integral.mapplot import mapplot
tmask = TheMask.mask_at_level(0)
fig, ax = mapplot({'data':tmask, 'clim':[0,1]})

ax.plot(I_red, J_red,'w.')
fig.show()









   



