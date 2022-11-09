import numpy as np
import scipy.io.netcdf as NC3
import netCDF4 as NC

mask3DVar = '/g100_scratch/userexternal/camadio0/RUN_ChlNitOxy_SUM/wrkdir/MODEL/DA_static_data/3D_VAR/GRID/BFM_gridFloat.nc'


MM = NC3.netcdf_file(mask3DVar,'r')
lon = MM.variables['lon'].data.copy()
lat = MM.variables['lat'].data.copy()
dep = MM.variables['dep'].data.copy()
lsm = MM.variables['lsm'].data.copy()
dx = MM.variables['dx'].data.copy()
dy = MM.variables['dy'].data.copy()
dz = MM.variables['dz'].data.copy()
regs = MM.variables['regs'].data.copy()
mdt = MM.variables['mdt'].data.copy()
tmsk = MM.variables['tmsk'].data.copy()

MM.close()

print('read mesh old')

dx_new = np.zeros_like(dx,dtype=np.float)
dy_new = np.zeros_like(dy,dtype=np.float)

longitude = lon.copy()
latitude = lat.copy()

radius = 6371.e3

km,jm,im = tmsk.shape


for j in range(jm):
    print(j)
    for i in range(1,im-1):
        dx_new[j,i] = radius*abs(longitude[j,i+1]-longitude[j,i-1])* \
                  0.5*np.pi/180. *np.cos(latitude[j,i]*np.pi/180.)
    dx_new[:,0] = dx_new[:,1]
    dx_new[:,im-1] = dx_new[:,im-2]

    for j in range(1,jm-1):
        for i in range(im):
            dy_new[j,i] = radius*abs(latitude[j+1,i]-latitude[j-1,i])* \
                      0.5*np.pi/180
    dy_new[0,:] = dy_new[1,:]
    dy_new[jm-1,:] = dy_new[jm-2,:]

print('new dx dy')

title='Mediterranean set-up of OPA model: MFS24125'

outfile = 'grid_dxdy_corrected.nc'
ncOUT = NC.Dataset(outfile,'w')

ncOUT.createDimension("im", im)
ncOUT.createDimension("jm", jm)
ncOUT.createDimension("km", km)

ncvar = ncOUT.createVariable('lon','f', ('jm','im'))
setattr(ncvar, 'units'        ,'degrees_east')
setattr(ncvar,'long_name'    ,'longitude')
setattr(ncvar, 'standard_name','longitude')
ncvar[:] = longitude

ncvar = ncOUT.createVariable('lat','f', ('jm','im'))
setattr(ncvar, 'units'        ,'degrees_north')
setattr(ncvar,'long_name'    ,'latitude')
setattr(ncvar, 'standard_name','latitude')
ncvar[:] = latitude

ncvar = ncOUT.createVariable('dep'   ,'f', ('km',))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'depth')
setattr(ncvar,'standard_name','depth')
ncvar[:] = dep

ncvar = ncOUT.createVariable('lsm'   ,'f', ('jm','im'))
setattr(ncvar,'long_name'    ,'Land-sea mask')
ncvar[:] = lsm


ncvar = ncOUT.createVariable('dx'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in X')
ncvar[:] = dx_new

ncvar = ncOUT.createVariable('dy'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in Y')
ncvar[:] = dy_new

ncvar = ncOUT.createVariable('dz'   ,'f', ('km',))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in Z')
ncvar[:] = dz

#ncvar = ncOUT.createVariable('topo'   ,'f', ('jm','im'))
#setattr(ncvar,'units'        ,'meters')
#setattr(ncvar,'long_name'    ,'Topography')
#ncvar[:] = topo

ncvar = ncOUT.createVariable('regs'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'')
setattr(ncvar,'long_name'    ,'Regions')
ncvar[:] = regs

ncvar = ncOUT.createVariable('mdt'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Mean dynamic topography')
ncvar[:] = mdt

ncvar = ncOUT.createVariable('tmsk'   ,'f', ('km','jm','im'))
setattr(ncvar,'units'        ,'')
setattr(ncvar,'long_name'    ,'T points mask')
ncvar[:] = tmsk

setattr(ncOUT,'title',title)

ncOUT.close()
