import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Create grid mask for Data assimilation
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output dir'''
                            )
    parser.add_argument(   '--name', '-n',
                            type = str,
                            required =False,
                            default = "DA_grid",
                            help = '''name of the output file'''
                            )
    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required =True,
                            default = "mask",
                            help = '''maskfile'''
                            )

    parser.add_argument(   '--case', '-c',
                            type = int,
                            required =False,
                            default = "3",
                            help = '''Define the calculation of the grid dimension: 1 = old way, 2 = 5000 m ,3 = today meshamask'''
                            )
    parser.add_argument(   '--basin', '-b',
                            type = str,
                            required =False,
                            default = "none",
                            help = '''if type subs the calculation will be on each subabsin defined, else on the whole mediterranean'''
                            )
    parser.add_argument(   '--gib', '-g',
                            type = str,
                            required =True,
                            help = '''Gibraltar index'''
                            )


    return parser.parse_args()

args = argument()


from bitsea.commons.mask import Mask
import numpy as np
from bitsea.commons import netcdf4
import netCDF4 as NC
from bitsea.commons.utils import addsep
from bitsea.commons.submask import SubMask
from bitsea.basins import V2



#mask=Mask("/gpfs/scratch/userexternal/ateruzzi/MASKS24_NRTV7C/meshmask.nc")
mask=Mask(args.maskfile)
CASE= args.case
name_file=args.name
outputdir=addsep(args.outdir)
basin_type=args.basin
LIMITE_OSCURAMENTO_GIB=float(args.gib)

km,jm,im =mask.shape

regs=np.zeros((jm,im),dtype=np.float32)
mdt=np.zeros((jm,im),dtype=np.float32)
longitude=mask.xlevels
latitude=mask.ylevels
depth=mask.zlevels
dx=np.zeros((jm,im),dtype=np.float32)
dy=np.zeros((jm,im),dtype=np.float32)
radius = 6371.e3

#CASE=2

if CASE == 1 :
#CASE 1 #OLD
    for j in range(jm):
        for i in range(1,im-1):
            dx[j,i]=radius*abs(longitude[j,i+1]-longitude[j,i-1])*0.5*np.pi/180. *np.cos(latitude[j,i]*np.pi/180.)
    dx[:,0]=dx[:,1]
    dx[:,im-1]=dx[:,im-2]

    for j in range(1,jm-1):
        for i in range(im):
            dy[j,i]=radius*abs(latitude[j+1,i]-latitude[j-1,i])*0.5*np.pi/180
    dy[0,:]=dy[1,:]
    dy[jm-1,:]=dy[jm-2,:]

elif CASE == 2 :

#CASE 2 #CASE 5000
    dx=np.ones((jm,im),dtype=np.float32)*5000
    dy=np.ones((jm,im),dtype=np.float32)*5000

else :
#CASE 3 #CASE mesh nostra
    dx=mask.e1t
    dy=mask.e2t


dz=mask.e3t[:,0,0]
#topo=mask.rough_bathymetry() #considering entire cell
#LIMITE_OSCURAMENTO_GIB=-4.9

tmsk=mask.mask.astype(np.float32)
for k in range(km):
    for j in range(jm):
        for i in range(im):
            if longitude[j,i] <= LIMITE_OSCURAMENTO_GIB:
                 tmsk[k,j,i]=0.0

lsm=tmsk[0,:,:]

basin_list=[V2.alb,V2.swm1,V2.swm2,V2.nwm,V2.tyr1,V2.tyr2,V2.adr,V2.aeg,V2.ion1,V2.ion2,V2.ion3,V2.lev1,V2.lev2,V2.lev3,V2.lev4]

mask0=mask.cut_at_level(0)
already_assigned=np.zeros(mask0.shape,dtype=np.bool)

if basin_type == 'subs':
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
    for j in range(jm):
        for i in range(im):
            if regs[j,i] == 0:
                regs[j,i] = 1        

else:
    k = 0
    for j in range(jm):
        for i in range(im):     
            if tmsk[1,j,i]:
                k = k + 1
                regs[j,i] = np.float32(k)
            else:
                regs[j,i] = 1.0


title='Mediterranean set-up of OPA model: MFS24125'

#outputdir = '/galileo/home/userexternal/gcoidess/CODICE_ANNA/preproc/DA/'
#name_file='prova'

outfile=outputdir + name_file+'.nc'
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

