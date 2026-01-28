import create_meshmask_nc as c_mask 

infile  = '/g100_work/OGS_test2528/camadio/RUN_RIVER_INLETS/MASK/ORIG_CMCC/mesh_mask.nc'

outfile = 'meshmask_CMCC.nc'

M       = c_mask.OrigMask(infile)
lon_cut = 0
depth_cut=0

Biscay_land = False
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land)

outfile = 'meshmask.nc'

lon_cut = M.getCutLocation(-8.875)
depth_cut=16
Biscay_land = True
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land)

M.nc_handler.close()
print('DONE!')
