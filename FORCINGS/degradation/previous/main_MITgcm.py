# SECTION FOR REDUCTION OF 1/16 meshmask
# nemo3.4 class reader
import create_meshmask_nc_MITgcm as c_mask

infile  = '/home/plazzari/data/ENEA/grid.t001.nc'

outfile = 'meshmask_MITgcm.nc'

M       = c_mask.OrigMask(infile)

free_surface=True
c_mask.create_meshmask_nc(M,outfile,free_surface)


M.nc_handler.close()
print('DONE!')
