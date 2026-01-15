# SECTION FOR REDUCTION OF 1/16 meshmask
# nemo3.4 class reader
import create_meshmask_nc_nemo as c_mask

infile  = '/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V4/mesh_mask_V1INGV.nc'

outfile = 'meshmask_INGVfor_ogstm.nc'

M       = c_mask.OrigMask(infile)

lon_cut = 0
depth_cut=0

free_surface=True
Biscay_land = False
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land,free_surface)

outfile = 'meshmask.nc'
 
lon_cut = 149 # or M.getCutLocation(-8.875)
depth_cut=2

Biscay_land = True
c_mask.create_meshmask_nc(M,outfile,lon_cut,depth_cut,Biscay_land,free_surface)


M.nc_handler.close()
print('DONE!')
