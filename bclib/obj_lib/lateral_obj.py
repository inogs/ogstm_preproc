import numpy as np
import numpy.matlib as npmat
#from scipy.io import netcdf as nc
import netCDF4 as nc
from bclib.io_lib import excel_obj as xlsobj
import logging
import code
import matplotlib.pyplot as plt

class lateral_bc:

    def __init__(self,mesh,file_nutrients):
        self.path = file_nutrients
        self._mesh_father = mesh
        self._extract_information()
        self._convert_information()
        self.season = (["0215-12:00:00","0515-12:00:00",
                        "0815-12:00:00","1115-12:00:00"])
        logging.info("lateral_bc builded")

    def _extract_information(self):
        try:
            self.ncfile = nc.Dataset(self.path, 'r')
        except:
            print("NUTRIENTS NOT FOUND")
            exit()

        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()

    def _convert_information(self):
        jpt = 4
        size_nutrients = np.zeros(4,dtype = np.int)
        size_nutrients[0] = jpt
        size_nutrients[1:] = self._mesh_father.tmask_dimension[1:]
        self.phos = np.zeros(size_nutrients)
        self.ntra = np.zeros(size_nutrients)
        self.dox = np.zeros(size_nutrients)
        self.sica = np.zeros(size_nutrients)
        self.dic = np.zeros(size_nutrients)
        self.alk = np.zeros(size_nutrients)
        tmask4d = np.zeros(size_nutrients)

        n_lev = self._mesh_father.tmask_dimension[1]
        vp_phos = np.zeros((4,n_lev));
        vp_ntra = np.zeros((4,n_lev));
        vp_sica = np.zeros((4,n_lev));
        vp_dox  = np.zeros((4,n_lev));
        vp_dic  = np.zeros((4,n_lev));
        vp_alk  = np.zeros((4,n_lev));

        nav_lev_in=np.concatenate(([0],self.lev1,[5500]))

        vp_phos_in=np.column_stack((self.N1p[:,1],self.N1p,self.N1p[:,len(self.lev1)-2]))

        vp_ntra_in=np.column_stack((self.N3n[:,1],self.N3n,self.N3n[:,len(self.lev1)-2]))

        vp_sica_in=np.column_stack((self.N5s[:,1],self.N5s,self.N5s[:,len(self.lev1)-2]))

        vp_dox_in=np.column_stack((self.O2o[ :,1],self.O2o,self.O2o[:,len(self.lev1)-2]))
        end = len(self.ALK)

        nav_lev_in2 = np.append([0],self.lev2.T)
        nav_lev_in2 = np.append(nav_lev_in2,[5500])

        vp_dic_in = np.append(self.DIC[1],self.DIC[:].T)
        vp_dic_in = np.append(vp_dic_in,self.DIC[end-1])

        vp_alk_in = np.append(self.ALK[1],self.ALK[:].T)
        vp_alk_in = np.append(vp_alk_in,self.ALK[end-1]);

        for i in range(4):
            jj= (~ np.isnan(vp_phos_in[i,:])) | ( vp_phos_in[i,:] < 1e19 ) ;
            vp_phos[i,:] = np.interp(self._mesh_father.nav_lev, nav_lev_in[jj], vp_phos_in[i,jj])
            jj=~ np.isnan(vp_ntra_in[i,:]) | ( vp_ntra_in[i,:] < 1e19 );
            vp_ntra[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_ntra_in[i,jj]);
            jj=~ np.isnan(vp_dox_in[i,:]) | ( vp_dox_in[i,:] < 1e19 );
            vp_dox[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_dox_in[ i,jj]);
            jj=~ np.isnan(vp_sica_in[i,:]) | ( vp_sica_in[i,:] < 1e19);
            vp_sica[i,:] = np.interp(self._mesh_father.nav_lev,nav_lev_in[jj],vp_sica_in[i,jj]);
            jj=~ np.isnan(vp_dic_in[:]) | ( vp_dic_in[:] < 1e19);
            vp_dic[i,:]  = np.interp(self._mesh_father.nav_lev,nav_lev_in2[jj],vp_dic_in[jj]);
            jj=~ np.isnan(vp_alk_in[:]) | ( vp_alk_in[:] < 1e19) ;
            vp_alk[i,:]  = np.interp(self._mesh_father.nav_lev,nav_lev_in2[jj],vp_alk_in[jj]);

#
# %Loop on time seasonal
            for jt in range(jpt):
                for jk in range(n_lev):
                    self.phos[jt,jk][:,:] = vp_phos[jt,jk];
                    self.ntra[jt,jk,:,:] = vp_ntra[jt,jk];
                    self.dox[jt,jk,:,:] = vp_dox[jt,jk];
                    self.sica[jt,jk,:,:] = vp_sica[jt,jk];
                    self.dic[jt,jk,:,:]  = vp_dic[jt,jk];
                    self.alk[jt,jk,:,:]  = vp_alk[jt,jk];


            for jt in range(jpt):
                    tmask4d[jt,:,:,:]=self._mesh_father.tmask[:,:,:];


            self.phos[tmask4d == 0]=1.e+20;
            self.ntra[tmask4d == 0]=1.e+20;
            self.dox[tmask4d == 0]=1.e+20;
            self.sica[tmask4d == 0]=1.e+20;
            self.dic[tmask4d == 0]=1.e+20;
            self.alk[ tmask4d == 0]=1.e+20;

