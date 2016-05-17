import create_meshmask_nc as c_mask 

infile  = 'mesh_mask_L121_zcr44_hth221_V3.4.nc'

outfile = 'prova_full.nc'

c_mask.create_meshmask_nc(infile,outfile,0,True)

outfile = 'prova_cut.nc'

c_mask.create_meshmask_nc(infile,outfile,300,True)

print('DONE!')
