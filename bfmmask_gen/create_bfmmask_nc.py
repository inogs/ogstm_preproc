import re
import numpy as np
import netCDF4 as nc

from commons.mask import Mask



def create_bfmmask_nc(meshmask_file="meshmask.nc", boundary_namelist="boundaries.nml", output_file="bfmmask.nc"):
    
    # create meshmask object
    meshmask = Mask.from_file(meshmask_file)
    
    # get dimensions
    (jpk, jpj, jpi) = meshmask.shape
    bfmmask = meshmask.mask
    
    # parse fortran namelist to get boundary infos and set mask accordingly
    with open(boundary_namelist, 'r') as boundaries:
        
        for b_line in boundaries:
            
            if (b_line[0] == '!'): continue
            
            b_infos = b_line.split()
            
            if (len(b_infos) > 1):
                
                b_type = b_infos[1].replace(",", "")
                b_file = b_infos[2].replace(",", "")
                
                # remove sponge points for every sponge boundary
                if (b_type == "SPO"): # only western sponge boundary supported

                    print("INFO: inside sponge branch")
                    
                    with open(b_file, 'r') as boundary:
                        for d_line in boundary:
                            if (re.search("length", d_line)):
                                length = float(d_line.split()[2].replace("d", "e"))
                                sample_lat = meshmask.ylevels[0,0]
                                (idx, _) = meshmask.convert_lon_lat_to_indices(lon=length, lat=sample_lat)
                                bfmmask[..., 0:idx+3] = False # hard-coded 3 cells of safety bound after the sponge
                
                # remove open boundary points for every open boundary
                elif (b_type == "OPE"):

                    print("INFO: inside open branch")

                    b_filenames_list = b_infos[3].replace(",", "")

                    with open(b_file, 'r') as boundary:
                        for d_line in boundary:
                            if (re.search("vars\(1\)", d_line)):
                                first_var = d_line.split()[2].replace("\"", "")

                    with open(b_filenames_list, 'r') as nc_files:
                        _ = nc_files.readline() # skip first line
                        sample_file = nc_files.readline().replace("\'", "").replace("\n", "")
                        sample_data = np.array(nc.Dataset(sample_file).variables[first_var])
                        b_points = np.logical_and(sample_data < 1.0e19, sample_data > 1.0e-6)
                        bfmmask = np.logical_and(bfmmask, ~b_points)
    
    # write to netcdf file
    out = nc.Dataset(output_file, 'w')
    
    out.createDimension("x", jpi);
    out.createDimension("y", jpj);
    out.createDimension("z", jpk);
    
    ncvar = out.createVariable("bfmmask", 'b', ("z", "y", "x"))
    
    ncvar[:] = bfmmask
    
    out.close()



if __name__ == "__main__":
    create_bfmmask_nc()
