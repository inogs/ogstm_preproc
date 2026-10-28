import create_meshmask_nc as c_mask 

# mesh_mask maust have all the masks and no masked land domains on glamt, gphit, etc.
mesh_file  = '/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc'

outfile = '/leonardo_work/RIVRD_OGS/med24v1_domain_cfg.nc'

#M       = c_mask.OrigMask(infile)
#lon_cut = 0
#depth_cut=0
#
#Biscay_land = False
#c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land)
#
#outfile = 'meshmask.nc'

M = c_mask.OrigMask(mesh_file)

#lat_cut, lon_cut = M.getCutLocation(-8.875, 36.909)
#lat_cut, lon_cut = M.getCutLocation(-8.83, 36.909)
lat_cut, lon_cut = M.getCutLocation(-8.7, 36.909)
depth_cut=16
Biscay_land = True
c_mask.create_meshmask_nc(M, outfile, lon_cut, depth_cut, Biscay_land)

#M.nc_handler.close()
print('DONE!')
