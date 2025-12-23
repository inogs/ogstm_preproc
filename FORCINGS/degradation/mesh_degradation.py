import numpy as np
import xarray as xr
from argparse import ArgumentParser
import yaml
from degrade_mesh import load_mesh, degrade_mesh, cut_med
from commons import dump_netcdf

'''
generates a reduced horizontal resolution mesh from an original mesh (with Atlantic buffer)
vertical resolution is unchanged.
resolution is degraded by merging cells in (ndeg x ndeg) blocks.
functining: 
0. load original mesh, 
1. expand lat, lon to make the new size a multiple of ndeg, extrapolate values properly
2. compute cell surface areas and volume (for weighted means)
3. degrade resolution of original variables
4. dump new meshmask.nc
usage: python degrade_mesh.py -y degrade_mesh.yaml

yaml file example:
maskfile: '/path/to/meshmask_in.nc'
outfile: '/path/to/meshmask_out.nc'
outfile_med: '/path/to/meshmask_out_med.nc'
ndeg: 6    #nr of cells to join in i, j
thresh: 1  #water point per block to assign to waterpoint in destination grid

'''

def argument():
    parser = ArgumentParser()
    parser.add_argument('--yamlfile','-y',
                        type=str,
                        required=True,
                        help='file with all parameters')
    return parser.parse_args()

def load_parameters():
    yamlfile=argument().yamlfile
    # yamlfile = 'degrade_mesh.yaml'
    with open(yamlfile) as f:
        Params = yaml.load(f, Loader=yaml.Loader)
    return(Params)

Params = load_parameters()
maskfile = Params['maskfile']
outfile = Params['outfile']
outfile_med = Params['outfile_med']
# by how much cells to degrade (nr of cells to join in i, j)
ndeg = Params['ndeg']
# threshold of points to assign waterpoint to degraded mesh
thresh = Params['thresh'] #THE ONLY OPTION THAT ENSURES NO RIVERS END UP ON LAND
# load meshmask and expand dimensions
print('Loading and expanding mesh...')
M1 = load_mesh(maskfile, ndeg)
#
print('Degrading mesh...')
M2 = degrade_mesh(M1, thresh, ndeg)
#
print('Dumping degraded 141 levels mesh ...')
dump_netcdf(M2, outfile)
#
print('Dumping degraded 125 levels mesh for ogstm ...')
MMed = cut_med(M2)
dump_netcdf(MMed, outfile_med)