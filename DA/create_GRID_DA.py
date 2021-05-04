

from commons.mask import Mask
import numpy as np
from commons import netcdf4
import netCDF4 as NC
from commons.utils import addsep
from commons.submask import SubMask
from basins import V2
import math

mask=Mask("/gpfs/scratch/userexternal/ateruzzi/MASKS24_NRTV7C/meshmask.nc")

km,jm,im =mask.shape

regs=np.zeros((jm,im),dtype=np.float32)
mdt=np.zeros((jm,im),dtype=np.float32)
longitude=mask.xlevels
latitude=mask.ylevels
depth=mask.zlevels

radius = 6371.e3
pi = math.acos(-1.)

CASE=2

if CASE == 1 :
#CASE 1 #OLD
    for j in range(jm):
        for i in range(1,im-1):
            dx[i,j]=radius*abs(longitude[i+1,j]-longitude[i-1,j])*0.5*pi/180. *cos(latitude[i,j]*pi/180.)
    dx[1,:]=dx[2,:]
    dx[im,:]=dx[im-1,:]

    for j in range(1,jm-1):
        for i in range(im):
            dy[i,j]=radius*abs(latitude[i,j+1]-latitude[i,j-1])*0.5*pi/180
    dy[:,1]=dy[:,2]
    dy[:jm]=dy[:,jm-1]

elif CASE == 2 :

#CASE 2 #CASE 5000
    dx=np.ones((jm,im),dtype=np.float32)*5000
    dy=np.ones((jm,im),dtype=np.float32)*5000

else :
#CASE 3 #CASE mesh nostra
    dx=mask.e1t
    dy=mask.e2t


dz=mask.e3t[:,0,0]
topo=mask.rough_bathymetry() #considering entire cell
LIMITE_OSCURAMENTO_GIB=-4.9

tmsk=mask.mask.astype(np.float32)
for k in range(km):
    for j in range(jm):
        for i in range(im):
            if longitude[j,i] <= LIMITE_OSCURAMENTO_GIB:
                 tmsk[k,j,i]=0.0

lsm=tmsk[0,:,:]

S=SubMask(V2.adr,maskobject=mask)

basin_list=[V2.alb,V2.swm1,V2.swm2,V2.nwm,V2.tyr1,V2.tyr2,V2.adr,V2.aeg,V2.ion1,V2.ion2,V2.ion3,V2.lev1,V2.lev2,V2.lev3,V2.lev4]

mask0=mask.cut_at_level(0)
already_assigned=np.zeros(mask0.shape,dtype=np.bool)

for isub,sub in enumerate(basin_list):
    S=SubMask(sub,maskobject=mask0)
    S.mask[already_assigned] = False
    already_assigned[S.mask] = True
    print sub.name
    S_matrix=S.mask[0,:,:].astype(np.float32)
    S_matrix=S_matrix * (isub+1)
    for j in range(jm):
        for i in range(im):
            #if S_matrix[j,i] == isub+1 :
            if S.mask[0,j,i] : 
                regs[j,i]=S_matrix[j,i]


title='Mediterranean set-up of OPA model: MFS24125'

outputdir = '/galileo/home/userexternal/gcoidess/CODICE_ANNA/preproc/DA/'
var_name='prova'

outfile=outputdir + var_name+'.nc'
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
ncvar[:] = depth

ncvar = ncOUT.createVariable('lsm'   ,'f', ('jm','im'))
setattr(ncvar,'long_name'    ,'Land-sea mask')
ncvar[:] = lsm

ncvar = ncOUT.createVariable('dx'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in X')
ncvar[:] = dx

ncvar = ncOUT.createVariable('dy'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in Y')
ncvar[:] = dy

ncvar = ncOUT.createVariable('dz'   ,'f', ('km',))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Resolution in Z')
ncvar[:] = dz

ncvar = ncOUT.createVariable('topo'   ,'f', ('jm','im'))
setattr(ncvar,'units'        ,'meters')
setattr(ncvar,'long_name'    ,'Topography')
ncvar[:] = topo

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

