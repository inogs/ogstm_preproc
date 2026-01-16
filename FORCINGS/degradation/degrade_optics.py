import numpy as np
import xarray as xr
import netCDF4
from argparse import ArgumentParser
from glob import glob
from pathlib import Path
#
import degrade_mesh as dm
from commons import degrade_wrap, load_parameters, reshape_blocks
import regridding as rg
from bitsea.commons.mask import Mask

'''
degrades optics input files:
aero.*.nc 
climatm.*.nc
.*.nc
takes the files already gridded on the 1/24 Med grid,
and degrades resolution by ndeg.
this is a dirty fix, if original ERA5 outputs are unavailable
'''

def argument():
    parser = ArgumentParser()
    parser.add_argument('--yamlfile','-y',
                        type=str,
                        required=True,
                        help='file with all parameters')
    return parser.parse_args()

def get_Area(mesh_in):
    '''
    gets t-cell [horizontal] surface area for weighted mean
    avoids loading full mesh, cause that's expansive
    '''
    M = xr.open_dataset(mesh_in)
    M1 = {}
    M1["e1t"] = dm.xpnd_wrap(M["e1t"], 'interp', ndeg)
    M1["e2t"] = dm.xpnd_wrap(M["e2t"], 'interp', ndeg)
    M1 = xr.Dataset(M1)
    M.close()
    #A0 = M1['e1t'].values[:] * M1['e2t'].values[:]
    A0 = M1['e1t'] * M1['e2t']
    return A0

def get_flist(prfx, Params, YYYY=None):
    '''
    this does all files at once
    consider looping on years for atm*.nc files
    '''
    indir = Params[f'indir_{prfx}']
    if YYYY is None:
        indir = indir.replace('YYYY', '????').replace('MM', '??') #for atm* files
    else:
        indir = indir.replace('YYYY', str(YYYY)).replace('MM', '??') #for atm* files
    infile = Params[f'infile_{prfx}']
    flist = sorted(glob(indir+infile))
    itrbl = [(ff.split('.')[-2], Path(ff)) for ff in flist]
    return itrbl

# indir = '/leonardo_work/OGS_devC_0/OPTICS/'
def load_xpnd_opt(infile, ndeg=1):
    '''
    loads optics file and expands it in 'edge' mode
    Arguments:
        infile : path to restart file
        ndeg : by how much to degrade resolution 
    Returns:
    DI : xarray Dataset, with 3D variables already expanded
         in 'edge' mode
    '''
    DI = {}
    D = xr.open_dataset(infile)
    vlist = list(D.variables.keys())
    vlist.remove('lat')
    vlist.remove('lon')
    #D = D.rename({'lon':'x', 'lat':'y', 'depth':'z_a'})
    D = D.rename({'lon':'x', 'lat':'y'})
    if 'depth' in D.dims:
        D = D.rename({'depth':'z_a'})
    #DI['lat'] = dm.xpnd_wrap(D['lat'], 'interp', ndeg)
    #DI['lon'] = dm.xpnd_wrap(D['lon'], 'interp', ndeg)
    for vname in vlist:
        A = D[vname]
        Ax = dm.xpnd_wrap(A, 'edge', ndeg)
        DI[vname] = Ax
    D.close()
    DI = xr.Dataset(DI)
    return DI

def awmean(X, W):
    '''
    W (time: 1, z_a: 1, y: 66, y_b: 6, x: 182, x_b: 6)
    mean, weighted by surface area or volume(t-grid)
    no mask, cause atmospheric fields
    also works on 3d if W has depth, as in aero* files
    '''
    weighted_sum = (W * X).sum(dim=('y_b', 'x_b'), skipna=True)
    weight_sum = W.sum(dim=('y_b', 'x_b'), skipna=True)
    wmean = weighted_sum / weight_sum
    assert np.isnan(wmean).any().values == False, "NaN values found in weighted mean!"
    return wmean

def climatm_writer(outfile, Od, TheMask):
    jpk, jpj, jpi = TheMask.shape
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);

    ncvar = ncOUT.createVariable('lon' ,'d',('lat','lon')); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('lat' ,'d',('lat','lon')); ncvar[:] = TheMask.ylevels

    setattr(ncOUT['lon'], 'units', 'degrees')
    setattr(ncOUT['lon'], 'long_name', 'longitude')
    setattr(ncOUT['lat'], 'units', 'degrees')
    setattr(ncOUT['lat'], 'long_name', 'longitude')

    ncvar = ncOUT.createVariable('cdrem', 'd', ('lat','lon')) ; ncvar[:] = Od['cdrem'].values[:]
    ncvar = ncOUT.createVariable('cldtcm', 'd', ('lat','lon')); ncvar[:] = Od['cldtcm'].values[:]

    setattr(ncOUT['cdrem'], 'units', '[um]')
    setattr(ncOUT['cdrem'], 'long_name', 'cloud droplet effective radius')
    setattr(ncOUT['cldtcm'], 'units', 'TODO')
    setattr(ncOUT['cldtcm'], 'long_name', 'TODO')

    ncOUT.close()
    return

def aero_writer(outfile, Od, TheMask):
    jpk, jpj, jpi = TheMask.shape
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);
    ncOUT.createDimension('depth',jpk);
    #
    ncvar = ncOUT.createVariable('lon' ,'d',('lat','lon')); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('lat' ,'d',('lat','lon')); ncvar[:] = TheMask.ylevels
    ncvar = ncOUT.createVariable('depth' ,'d',('depth')); ncvar[:] = TheMask.zlevels
    #
    setattr(ncOUT['lon'], 'units', 'degrees')
    setattr(ncOUT['lon'], 'long_name', 'longitude')
    setattr(ncOUT['lat'], 'units', 'degrees')
    setattr(ncOUT['lat'], 'long_name', 'longitude')
    setattr(ncOUT['depth'], 'long_name', 'depth')
    setattr(ncOUT['depth'], 'units', '[m]')
    #
    vlist = list(Od.variables.keys())
    for var in vlist:
        ncvar = ncOUT.createVariable(var, 'd', ('depth','lat','lon')) ; ncvar[:] = Od[var].values[:]
    #
    setattr(ncOUT['taua'], 'long_name', 'aerosol optical thickness')
    setattr(ncOUT['taua'], 'units', '[-]')
    setattr(ncOUT['taua'], 'orig', 'MODIS_AEROSOL')
    setattr(ncOUT['asymp'], 'long_name', 'aerosol asymmetry parameter')
    setattr(ncOUT['asymp'], 'units', '[-]')
    setattr(ncOUT['asymp'], 'orig', 'MODIS_AEROSOL')
    setattr(ncOUT['ssalb'], 'long_name', 'aerosol single scattering albedo')
    setattr(ncOUT['ssalb'], 'units', '[-]')
    setattr(ncOUT['ssalb'], 'orig', 'MODIS_AEROSOL')
    #
    ncOUT.close()
    return

def atm_writer(outfile, Od, TheMask):
    jpk, jpj, jpi = TheMask.shape
    ncOUT=netCDF4.Dataset(outfile,"w", format="NETCDF4")
    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);
    #
    ncvar = ncOUT.createVariable('lon' ,'d',('lat','lon')); ncvar[:] = TheMask.xlevels
    ncvar = ncOUT.createVariable('lat' ,'d',('lat','lon')); ncvar[:] = TheMask.ylevels
    #
    vlist = list(Od.variables.keys())
    for var in vlist:
        ncvar = ncOUT.createVariable(var, 'd', ('lat','lon')) ; ncvar[:] = Od[var].values[:]
    #
    setattr(ncOUT['lon'], 'units', 'degrees')
    setattr(ncOUT['lon'], 'long_name', 'longitude')    
    setattr(ncOUT['lat'], 'units', 'degrees')
    setattr(ncOUT['lat'], 'long_name', 'latitude')
    setattr(ncOUT['sp'], 'long_name', 'surface pressure')
    setattr(ncOUT['sp'], 'units', 'Pa')
    setattr(ncOUT['sp'], 'orig', 'ERA5')
    setattr(ncOUT['msl'], 'long_name', 'mean sea level pressure')
    setattr(ncOUT['msl'], 'units', 'Pa')
    setattr(ncOUT['msl'], 'orig', 'ERA5')
    setattr(ncOUT['u10'], 'long_name', 'zonal wind velocity')
    setattr(ncOUT['u10'], 'units', 'm s*-1')
    setattr(ncOUT['u10'], 'orig', 'ERA5')
    setattr(ncOUT['v10'], 'long_name', 'meridional wind velocity')
    setattr(ncOUT['v10'], 'units', 'm s*-1')
    setattr(ncOUT['v10'], 'orig', 'ERA5')
    setattr(ncOUT['tclw'], 'long_name', 'Total column cloud liquid water')
    setattr(ncOUT['tclw'], 'units', 'kg m**-2')
    setattr(ncOUT['tclw'], 'orig', 'ERA5')
    setattr(ncOUT['tco3'], 'long_name', 'Total column ozone')
    setattr(ncOUT['tco3'], 'units', 'kg m**-2')
    setattr(ncOUT['tco3'], 'orig', 'ERA5')
    setattr(ncOUT['t2m'], 'long_name', '2 metre temperature')
    setattr(ncOUT['t2m'], 'units', 'K')
    setattr(ncOUT['t2m'], 'orig', 'ERA5')
    setattr(ncOUT['d2m'], 'long_name', '2 metre dewpoint temperature')
    setattr(ncOUT['d2m'], 'units', 'K')
    setattr(ncOUT['d2m'], 'orig', 'ERA5')
    setattr(ncOUT['tcc'], 'long_name', 'Total cloud cover')
    setattr(ncOUT['tcc'], 'units', '[-]')
    setattr(ncOUT['tcc'], 'orig', 'ERA5')
    #
    ncOUT.close()
    return

def degrade_opt(DI, A0, ndeg):
    '''
    NB, despite the naming also works on 3D fields!
    '''
    #DI = DI.rename({'lon':'x', 'lat':'y', 'depth':'z_a'})
    #DI = DI.rename({'lon':'x', 'lat':'y'})
    Od = {}
    vlist = list(DI.variables.keys())
    vlist.remove('y')
    vlist.remove('x')
    for var2d in vlist:
        X = reshape_blocks(DI[var2d], ndeg)
        Od[var2d] = awmean(X, A0)
    Od = xr.Dataset(Od)
    Od = Od.rename({'x':'lon', 'y':'lat'})
    if 'z' in Od.dims:
        Od = Od.rename({'z_a':'depth'})
    Od = Od.squeeze() #remove dims 'time' and 'z_a'
    return Od 

def make_outdir(outfile:Path):
    '''
    Generates /yyyy/mm/ subdirectory in outdir based on outfile name
    and creates it if it does not exist
    '''
    outdir = outfile.parent
    outdir.mkdir(exist_ok=True, parents=True)
    return

if __name__=='__main__':
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nranks = comm.size
    except:
        comm = None
        rank = 0
        nranks = 1
    # yamlfile = 'degrade_optics.yaml'
    Params = load_parameters(argument().yamlfile)
    mesh_in = Path(Params['mesh_in'])
    mesh_out = Path(Params['mesh_out'])
    outdir = Path(Params['outdir'])
    ndeg = Params['ndeg']
    atm_y0 = Params['atm_y0']
    atm_yE = Params['atm_yE']

    do_atm = Params['do_atm']
    do_climatm = Params['do_climatm']
    do_aero = Params['do_aero']

    A0 = dm.reshape_blocks(get_Area(mesh_in), ndeg)
    Maskout = Mask.from_file(mesh_out)
    A03d = (A0.isel(z_a=0).expand_dims(z_a=Maskout.shape[0])).transpose(*A0.dims)

    itrbl_aero = get_flist('aero', Params)
    itrbl_climatm = get_flist('climatm', Params)

    if do_aero:
        print('--- DEGRADING aero FILES ---')
        for dt_str, fname in itrbl_aero[rank::nranks]:
            print(fname.name)
            outfile = outdir / fname.name
            make_outdir(outfile)
            DI = load_xpnd_opt(fname, ndeg)
            Od = degrade_opt(DI, A03d, ndeg)
            aero_writer(outfile, Od, Maskout)

        # Wait until all ranks reach this point (just to have the output clean)
        if nranks>1: 
            comm.Barrier()

    if do_climatm:
        print('--- DEGRADING climatm FILES ---')
        for dt_str, fname in itrbl_climatm[rank::nranks]:
            print(fname.name)
            outfile = outdir / fname.name
            make_outdir(outfile)
            DI = load_xpnd_opt(fname, ndeg)
            Od = degrade_opt(DI, A0, ndeg)
            climatm_writer(outfile, Od, Maskout)

        if nranks>1: 
            comm.Barrier()

    if do_atm:
        print('--- DEGRADING atm FILES ---')
        for yyyy in range(atm_y0, atm_yE+1):
            itrbl_atm = get_flist('atm', Params, yyyy)
            for dt_str, fname in itrbl_atm[rank::nranks]:
                print(fname.name)
                outfile = outdir / fname.parts[-3] / fname.parts[-2] / fname.name
                make_outdir(outfile)
                DI = load_xpnd_opt(fname, ndeg)
                Od = degrade_opt(DI, A0, ndeg)
                atm_writer(outfile, Od, Maskout)

