
import numpy as np
import scipy.io.netcdf as nc


class boun_mesh:

    """
        Class bounmesh

    """

    def __init__(self, ncfile):
        self.path = ncfile




class sub_mesh:
    """
        Class sub mesh

    """

    def __init__(self, ncfile):
        self.path = ncfile
        self.extract_information()


    def extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            setattr(self, i, self.ncfile.variables[i])


class mesh:
    """
        Class mesh

    """

    def __init__(self, ncf_mesh,ncf_submesh,ncf_bounmesh):
        self.path = ncf_mesh
        self.extract_information()
        self.submesh = sub_mesh(ncf_submesh)
        self.bounmesh = boun_mesh(ncf_bounmesh)

    def extract_information(self):
        self.ncfile = nc.netcdf_file(self.path, 'r')
        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            setattr(self, i, self.ncfile.variables[i])

    def generate_boundmask(self,elab):

        """

        This fuction generate boundmask
        from elab we use :
        vnudg  : array of variable
        end_nudging : elab.end_nudging
        rdmpmin : elab.rdmpmin
        rdmpmax : elab.rdmpmin

        """

        vnudg = elab.varibles
        rdpmin = elab.rdpmin
        rdpmax = elab.rdpmax
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
