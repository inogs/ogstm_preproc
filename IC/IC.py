import netCDF4
def RSTwriter(outfile, var,rst,MaskObj):
    jpk, jpj, jpi = MaskObj.shape
    rst[~MaskObj.mask] = 1.e+20
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',1)

    TRN   = 'TRN' + var;
    ncvar = ncOUT.createVariable('nav_lon' ,'d',('y','x')           ); ncvar[:] = MaskObj.xlevels
    ncvar = ncOUT.createVariable('nav_lat' ,'d',('y','x')           ); ncvar[:] = MaskObj.ylevels
    ncvar = ncOUT.createVariable('nav_lev' ,'d',('z')               ); ncvar[:] = MaskObj.zlevels
    ncvar = ncOUT.createVariable('time'    ,'d',('time',)           ); ncvar    = 1.;
    ncvar = ncOUT.createVariable(TRN       ,'d',('time','z','y','x'),zlib=True); ncvar[:] = rst;


    setattr(ncOUT.variables[TRN]   ,'missing_value',1e+20                              );
    setattr(ncOUT.variables['time'],'Units'        ,'seconds since 1582-10-15 00:00:00');
    setattr(ncOUT                  ,'TimeString'   ,'20010101-00:00:00');
    ncOUT.close()