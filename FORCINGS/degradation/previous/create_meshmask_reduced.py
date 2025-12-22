import reducer
import numpy as np
isFreeSurface=True

inputmesh="meshmask.nc"

R =reducer.reducer(isFreeSurface)
R.read_fine(inputmesh)
passo=4
Mred,I,J = R.reduce(passo)
Mexp,I,J = R.expand(Mred,I,J)


outputmesh="meshmask_472.nc"
R.dumpfile(Mexp, outputmesh)

np.savetxt('South_West_I_indexes.txt',I,'%d') # of mesh expanded
np.savetxt('South_West_J_indexes.txt',J,'%d')

#Mexp= reducer.mesh("meshmask_472.nc")
M_med, cut = R.cut_med(Mexp, lon_cut=-8.875, depth_cut=2,biscay_land=True)
R.dumpfile(M_med, "meshmask_470.nc")

# For BC - input for river_reindexer
np.savetxt('ogstm_South_West_I_indexes.txt',I[cut:],'%d') # of M_med mesh
np.savetxt('ogstm_South_West_J_indexes.txt',J,'%d')
