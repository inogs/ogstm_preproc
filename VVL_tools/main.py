import calculate_scale_factors as c_mask 

infile  = 'mesh_mask.nc'

outfile = 'prova_mesh.nc'

M       = c_mask.OrigMask(infile)

lon_cut = M.getCutLocation(-8.875)
depth_cut=11

Biscay_land = True
free_surface=True

c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land, free_surface)

M.nc_handler.close()
print('DONE!')
