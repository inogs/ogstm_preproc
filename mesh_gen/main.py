import create_meshmask_nc as c_mask 

infile  = '/pico/scratch/userexternal/plazzari/TEST_24/mesh_mask_L121_zcr44_hth221_V3.4.nc'

outfile = 'prova_full.nc'

M       = c_mask.OrigMask(infile)
lon_cut = 0
depth_cut=0

free_surface=True

c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,free_surface)

outfile = 'prova_cut.nc'

lon_cut = M.getCutLocation(-8.875)
depth_cut=5

c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,free_surface)

M.nc_handler.close()
print('DONE!')
