import numpy as np
import netCDF4 as nc
import logging



class lateral_bc:

    def __init__(self,conf,maskobj):
        self.path = conf.file_nutrients
        self._extract_information()
        self.season = conf.gib_season
        self._convert_information(maskobj)
        logging.info("lateral_bc builded")

    def _extract_information(self):
        self.ncfile = nc.Dataset(self.path, 'r')

        for i in self.ncfile.dimensions:
            setattr(self, i, self.ncfile.dimensions[i])
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()

    def _convert_information(self,mask):
        jpt = len(self.season)
        print('self.season',self.season)
        print('--------------------')
        jpk, jpj, jpi = mask.shape
        size_nutrients = (jpt,jpk,jpj, jpi)
        self.phos = np.zeros(size_nutrients, np.float64)
        self.ntra = np.zeros(size_nutrients, np.float64)
        self.dox  = np.zeros(size_nutrients, np.float64)
        self.sica = np.zeros(size_nutrients, np.float64)
        self.dic  = np.zeros(size_nutrients, np.float64)
        self.alk  = np.zeros(size_nutrients, np.float64)
        self.me0  = np.zeros(size_nutrients, np.float64)
        self.me2  = np.zeros(size_nutrients, np.float64)
        self.mmh  = np.zeros(size_nutrients, np.float64)
        self.dmh  = np.zeros(size_nutrients, np.float64)
        self.pl1  = np.zeros(size_nutrients, np.float64)
        self.pl2  = np.zeros(size_nutrients, np.float64)
        self.pl3  = np.zeros(size_nutrients, np.float64)
        self.pl4  = np.zeros(size_nutrients, np.float64)
        self.zo6  = np.zeros(size_nutrients, np.float64)
        self.zo5  = np.zeros(size_nutrients, np.float64)
        self.zo4  = np.zeros(size_nutrients, np.float64)
        self.zo3  = np.zeros(size_nutrients, np.float64)
        tmask4d   = np.zeros(size_nutrients, np.bool_)
      # tmask4d   = np.zeros(size_nutrients, np.bool)
        nav_lev = mask.zlevels
        n_lev = jpk
        print('n_lev',n_lev)
        #n_lev2 = jpk
        vp_phos = np.zeros((jpt,n_lev));
        vp_ntra = np.zeros((jpt,n_lev));
        vp_sica = np.zeros((jpt,n_lev));
        vp_dox  = np.zeros((jpt,n_lev));
        vp_dic  = np.zeros((jpt,n_lev));
        vp_alk  = np.zeros((jpt,n_lev));
        vp_Hg0  = np.zeros((jpt,n_lev));
        vp_Hg2  = np.zeros((jpt,n_lev));
        vp_MHg = np.zeros((jpt,n_lev));
        vp_DHg = np.zeros((jpt,n_lev));
        vp_P1h = np.zeros((jpt,n_lev));
        vp_P2h = np.zeros((jpt,n_lev));
        vp_P3h = np.zeros((jpt,n_lev));
        vp_P4h = np.zeros((jpt,n_lev));
        vp_Z6h = np.zeros((jpt,n_lev));
        vp_Z5h = np.zeros((jpt,n_lev));
        vp_Z4h = np.zeros((jpt,n_lev));
        vp_Z3h = np.zeros((jpt,n_lev));

        nav_lev_in=np.concatenate(([0],self.lev1,[5500]))

        vp_phos_in=np.column_stack((self.N1p[:,0], self.N1p, self.N1p[:,-2]))
        vp_ntra_in=np.column_stack((self.N3n[:,0], self.N3n, self.N3n[:,-2]))
        vp_sica_in=np.column_stack((self.N5s[:,0], self.N5s, self.N5s[:,-2]))
        vp_dox_in =np.column_stack((self.O2o[:,0], self.O2o, self.O2o[:,-2]))
#######  vp_Hg0_in =np.column_stack((self.Hg0[:,0], self.Hg0, self.Hg0[:,-2]))
#### end of two dimensional vars - 1D vars

        end = len(self.ALK)
        print('END',end)
        nav_lev_in2 = np.append([0],self.lev2.T)
        nav_lev_in2 = np.append(nav_lev_in2,[5500])
        print('NAVin2',nav_lev_in2 ) 
        vp_dic_in = np.append(self.DIC[0],self.DIC[:].T)
        vp_dic_in = np.append(vp_dic_in,self.DIC[end-1])

        vp_alk_in = np.append(self.ALK[0],self.ALK[:].T)
        vp_alk_in = np.append(vp_alk_in,self.ALK[end-1]);

        print('vp_alk_in0',vp_alk_in)

        end2=len(self.Hg0)
        print('end for Hg0',end2)
        nav_lev_in3 = self.lev3.T
        print('NAVin3',nav_lev_in3) 

        vp_Hg0_in = self.Hg0[:].T
       # vp_Hg0  = np.zeros((jpt,n_lev))
        print('self.Hg0',self.Hg0)
     #   print('vp_Hg0_in1',vp_Hg0_in)
     #   print('self.Hg0[0]',self.Hg0[0])
     #   print('self.Hg0[0].T',self.Hg0[0].T)

        vp_Hg2_in = self.Hg2[:].T
        vp_MHg_in = self.MMHg[:].T
        vp_DHg_in = self.DMHg[:].T

        for i in range(jpt):
            jj= (~ np.isnan(vp_phos_in[i,:])) | ( vp_phos_in[i,:] < 1e19 )
            vp_phos[i,:] = np.interp(nav_lev, nav_lev_in[jj], vp_phos_in[i,jj])
            jj=~ np.isnan(vp_ntra_in[i,:]) | ( vp_ntra_in[i,:] < 1e19 )
            vp_ntra[i,:] = np.interp(nav_lev,nav_lev_in[jj],vp_ntra_in[i,jj])
            jj=~ np.isnan(vp_dox_in[i,:]) | ( vp_dox_in[i,:] < 1e19 )
            vp_dox[i,:] = np.interp(nav_lev,nav_lev_in[jj],vp_dox_in[ i,jj])
            jj=~ np.isnan(vp_sica_in[i,:]) | ( vp_sica_in[i,:] < 1e19)
            vp_sica[i,:] = np.interp(nav_lev,nav_lev_in[jj],vp_sica_in[i,jj])
            jj=~ np.isnan(vp_dic_in[:]) | ( vp_dic_in[:] < 1e19)
            vp_dic[i,:]  = np.interp(nav_lev,nav_lev_in2[jj],vp_dic_in[jj])
            jj=~ np.isnan(vp_alk_in[:]) | ( vp_alk_in[:] < 1e19)
            vp_alk[i,:]  = np.interp(nav_lev,nav_lev_in2[jj],vp_alk_in[jj])

            jj=~ np.isnan(vp_Hg0_in[:]) | ( vp_Hg0_in[:] < 1e19)
            vp_Hg0[i,:]  = np.interp(nav_lev,nav_lev_in3[jj],vp_Hg0_in[jj])
            jj=~ np.isnan(vp_Hg2_in[:]) | ( vp_Hg2_in[:] < 1e19)
            vp_Hg2[i,:]  = np.interp(nav_lev,nav_lev_in3[jj],vp_Hg2_in[jj])
            jj=~ np.isnan(vp_MHg_in[:]) | ( vp_MHg_in[:] < 1e19)
            vp_MHg[i,:]  = np.interp(nav_lev,nav_lev_in3[jj],vp_MHg_in[jj])
            jj=~ np.isnan(vp_DHg_in[:]) | ( vp_DHg_in[:] < 1e19)
            vp_DHg[i,:]  = np.interp(nav_lev,nav_lev_in3[jj],vp_DHg_in[jj])

            print('vp_Hg0_in[:]',vp_Hg0_in[:].shape)
            print('vp_alk_in[:]',vp_alk_in[:].shape)

            print('self.phos',self.phos.shape)
            print('self.phos',self.phos.shape)
            print('jpts',jpt)
            print('n_levs',n_lev)

            print('vp.phos',vp_phos.shape)
            print('vp_Hg0',vp_Hg0.shape)
            print('tmask4d',tmask4d.shape)
# %Loop on time seasonal
            for jt in range(jpt):
                for jk in range(n_lev):
                    self.phos[jt,jk,:,:] = vp_phos[jt,jk]
                    self.ntra[jt,jk,:,:] = vp_ntra[jt,jk]
                    self.dox[ jt,jk,:,:] = vp_dox[ jt,jk]
                    self.sica[jt,jk,:,:] = vp_sica[jt,jk]
                    self.dic[ jt,jk,:,:] = vp_dic[ jt,jk]
                    self.alk[ jt,jk,:,:] = vp_alk[ jt,jk]
                    self.me0[ jt,jk,:,:] = vp_Hg0[ jt,jk]
                    self.me2[ jt,jk,:,:] = vp_Hg2[ jt,jk]
                    self.mmh[ jt,jk,:,:] = vp_MHg[ jt,jk]
                    self.dmh[ jt,jk,:,:] = vp_DHg[ jt,jk]
                    self.pl1[ jt,jk,:,:] = vp_P1h[ jt,jk]
                    self.pl2[ jt,jk,:,:] = vp_P2h[ jt,jk]
                    self.pl3[ jt,jk,:,:] = vp_P3h[ jt,jk]
                    self.pl4[ jt,jk,:,:] = vp_P4h[ jt,jk]
                    self.zo6[ jt,jk,:,:] = vp_Z6h[ jt,jk]
                    self.zo5[ jt,jk,:,:] = vp_Z5h[ jt,jk]
                    self.zo4[ jt,jk,:,:] = vp_Z4h[ jt,jk]
                    self.zo3[ jt,jk,:,:] = vp_Z3h[ jt,jk]


            for jt in range(jpt): tmask4d[jt,:,:,:]=mask[:]

            self.phos[~tmask4d]=1.e+20
            self.ntra[~tmask4d]=1.e+20
            self.dox[ ~tmask4d]=1.e+20
            self.sica[~tmask4d]=1.e+20
            self.dic[ ~tmask4d]=1.e+20
            self.alk[ ~tmask4d]=1.e+20
            self.me0[ ~tmask4d]=1.e+20
            self.me2[ ~tmask4d]=1.e+20
            self.mmh[ ~tmask4d]=1.e+20
            self.dmh[ ~tmask4d]=1.e+20
