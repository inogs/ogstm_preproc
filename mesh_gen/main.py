import create_meshmask_nc as c_mask 

infile  = '/pico/scratch/userexternal/plazzari/TEST_24/FORCINGS/STATIC_DATA/mesh_mask_L141_zcr64_hth111_r8_V3.6.nc'

outfile = 'prova_full.nc'

M       = c_mask.OrigMask(infile)
lon_cut = 0
depth_cut=0

free_surface=True
Biscay_land = False
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land,free_surface)

outfile = '/pico/scratch/userexternal/plazzari/TEST_24/FORCINGS/STATIC_DATA/meshmask_FS.nc'

lon_cut = M.getCutLocation(-8.875)
depth_cut=11
Biscay_land = True
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land, free_surface)

M.nc_handler.close()
print('DONE!')
