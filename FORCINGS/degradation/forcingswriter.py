import netCDF4 as NC
from bitsea.commons.mask import Mask
#Maskfile_OGS="meshmask_6125.nc"
#TheMask = Mask(Maskfile_OGS)

def writefileU(filename, mask, U, taux):
    
    jpk,jpj,jpi = mask.shape
    
    
    ncOUT = NC.Dataset(filename,"w",format="NETCDF4")
    ncOUT.createDimension('x',jpi)
    ncOUT.createDimension('y',jpj)
    ncOUT.createDimension('depthu',jpk)
    ncOUT.createDimension('time_counter',1)
    
    ncvar = ncOUT.createVariable('nav_lon','f',('y','x'))
    ncvar[:]=mask.xlevels
    ncvar = ncOUT.createVariable('nav_lat','f',('y','x'))
    ncvar[:]=mask.ylevels
    ncvar = ncOUT.createVariable('depthu','f',('depthu',))
    ncvar[:]=mask.zlevels
    
    ncvar=ncOUT.createVariable('vozocrtx', 'f' ,('time_counter', 'depthu', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:]=U
    setattr(ncvar,'standard_name',"sea_water_x_velocity")
    setattr(ncvar,'long_name',"sea_water_x_velocity")
    setattr(ncvar,'units',"m/s" )
    setattr(ncvar,'online_operation',"instant")
    setattr(ncvar,'interval_operation',"1 d")
    setattr(ncvar,'interval_write',"1 d")
    setattr(ncvar,'cell_methods',"time: point")
    setattr(ncvar,'coordinates',"time_instant depthu nav_lon nav_lat")
    
    ncvar = ncOUT.createVariable('sozotaux', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:]=taux
    setattr(ncvar,'standard_name',"surface_downward_x_stress")
    setattr(ncvar,'long_name',"surface_downward_x_stress")
    setattr(ncvar,'units',"N/m2" )
    ncOUT.close()
    
def writefileV(filename, mask, V, tauy):
    
    jpk,jpj,jpi = mask.shape
    
    
    ncOUT = NC.Dataset(filename,"w",format="NETCDF4")
    ncOUT.createDimension('x',jpi)
    ncOUT.createDimension('y',jpj)
    ncOUT.createDimension('depthv',jpk)
    ncOUT.createDimension('time_counter',1)
    
    ncvar = ncOUT.createVariable('nav_lon','f',('y','x'))
    ncvar[:]=mask.xlevels
    ncvar = ncOUT.createVariable('nav_lat','f',('y','x'))
    ncvar[:]=mask.ylevels
    ncvar = ncOUT.createVariable('depthv','f',('depthv',))
    ncvar[:]=mask.zlevels
    
    ncvar=ncOUT.createVariable('vomecrty', 'f' ,('time_counter', 'depthv', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:]=V
    setattr(ncvar,'standard_name',"sea_water_y_velocity")
    setattr(ncvar,'long_name',"sea_water_y_velocity")
    setattr(ncvar,'units',"m/s" )
    setattr(ncvar,'online_operation',"instant")
    setattr(ncvar,'interval_operation',"1 d")
    setattr(ncvar,'interval_write',"1 d")
    setattr(ncvar,'cell_methods',"time: point")
    setattr(ncvar,'coordinates',"time_instant depthv nav_lon nav_lat")
    
    ncvar = ncOUT.createVariable('sometauy', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:]=tauy
    setattr(ncvar,'standard_name',"surface_downward_y_stress")
    setattr(ncvar,'long_name',"surface_downward_y_stress")
    setattr(ncvar,'units',"N/m2" )
    ncOUT.close()

def writefileW(filename, mask, W, K):
    
    jpk,jpj,jpi = mask.shape
    
    
    ncOUT = NC.Dataset(filename,"w",format="NETCDF4")
    ncOUT.createDimension('x',jpi)
    ncOUT.createDimension('y',jpj)
    ncOUT.createDimension('depthw',jpk)
    ncOUT.createDimension('time_counter',1)
    
    ncvar = ncOUT.createVariable('nav_lon','f',('y','x'))
    ncvar[:]=mask.xlevels
    ncvar = ncOUT.createVariable('nav_lat','f',('y','x'))
    ncvar[:]=mask.ylevels
    ncvar = ncOUT.createVariable('depthw','f',('depthw',))
    ncvar[:]=mask.zlevels
    if W is not None:
        ncvar=ncOUT.createVariable('vovecrtz', 'f' ,('time_counter', 'depthw', 'y', 'x'),fill_value=1.0e+20)
        ncvar[:]=W
        setattr(ncvar,'standard_name',"sea_water_z_velocity")
        setattr(ncvar,'long_name',"sea_water_z_velocity")
        setattr(ncvar,'units',"m/s" )
        setattr(ncvar,'online_operation',"instant")
        setattr(ncvar,'interval_operation',"1 d")
        setattr(ncvar,'interval_write',"1 d")
        setattr(ncvar,'cell_methods',"time: point")
        setattr(ncvar,'coordinates', "time_instant depthu nav_lon nav_lat")
    
    
    ncvar=ncOUT.createVariable('votkeavt', 'f' ,('time_counter', 'depthw', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:]=K
    
    setattr(ncvar,'standard_name',"sea_water_z_velocity")
    setattr(ncvar,'long_name',"sea_water_z_velocity")
    setattr(ncvar,'units',"m/s" )
    setattr(ncvar,'online_operation',"instant")
    setattr(ncvar,'interval_operation',"1 d")
    setattr(ncvar,'interval_write',"1 d")
    setattr(ncvar,'cell_methods',"time: point")

    setattr(ncvar,'coordinates',"time_instant depthu nav_lon nav_lat")
    

    ncOUT.close()

def writefileT(filename, mask, T, S, ETA, SW_Rad, RUNOFF, MXL):
    
    jpk,jpj,jpi = mask.shape
    
    
    ncOUT = NC.Dataset(filename,"w",format="NETCDF4")
    ncOUT.createDimension('x',jpi)
    ncOUT.createDimension('y',jpj)
    ncOUT.createDimension('deptht',jpk)
    ncOUT.createDimension('time_counter',1)
    
    ncvar = ncOUT.createVariable('nav_lon','f',('y','x'))
    ncvar[:] = mask.xlevels
    ncvar = ncOUT.createVariable('nav_lat','f',('y','x'))
    ncvar[:] = mask.ylevels
    ncvar = ncOUT.createVariable('deptht','f',('deptht',))
    ncvar[:] = mask.zlevels
    
    ncvar = ncOUT.createVariable('vosaline', 'f' ,('time_counter', 'deptht', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = S
    setattr(ncvar,'standard_name',"sea_water_z_velocity")
    setattr(ncvar,'long_name',"sea_water_z_velocity")
    setattr(ncvar,'units',"m/s" )
    setattr(ncvar,'online_operation',"instant")
    setattr(ncvar,'interval_operation',"1 d")
    setattr(ncvar,'interval_write',"1 d")
    setattr(ncvar,'cell_methods',"time: point")
    setattr(ncvar,'coordinates', "time_instant depthu nav_lon nav_lat")
    
    
    ncvar = ncOUT.createVariable('votemper', 'f' ,('time_counter', 'deptht', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = T
    
    setattr(ncvar,'standard_name',"sea_water_z_velocity")
    setattr(ncvar,'long_name',"sea_water_z_velocity")
    setattr(ncvar,'units',"m/s" )
    setattr(ncvar,'online_operation',"instant")
    setattr(ncvar,'interval_operation',"1 d")
    setattr(ncvar,'interval_write',"1 d")
    setattr(ncvar,'cell_methods',"time: point")
    setattr(ncvar,'coordinates',"time_instant depthu nav_lon nav_lat")
    
    ncvar = ncOUT.createVariable('sossheig', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = ETA
    setattr(ncvar,'standard_name',"sea_surface_height_above_geoid")
    setattr(ncvar,'long_name',"sea_surface_height_above_geoid")
    setattr(ncvar,'units',"m" )

    ncvar = ncOUT.createVariable('soshfldo', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = SW_Rad
    setattr(ncvar,'standard_name',"Shortwave Radiation")
    setattr(ncvar,'long_name',"Shortwave Radiation")
    setattr(ncvar,'units',"W/m2" )

    ncvar = ncOUT.createVariable('sorunoff', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = RUNOFF
    setattr(ncvar,'standard_name',"water_flux_into_sea_water_from_rivers")
    setattr(ncvar,'long_name',"water_flux_into_sea_water_from_rivers")
    setattr(ncvar,'units',"kg/m2/s" )

    ncvar = ncOUT.createVariable('somxl010', 'f', ('time_counter', 'y', 'x'),fill_value=1.0e+20)
    ncvar[:] = MXL
    setattr(ncvar,'standard_name',"mixed_layer_depth_0.01")
    setattr(ncvar,'long_name',"mixed_layer_depth_0.01")
    setattr(ncvar,'units',"m" )

    ncOUT.close()
   
