import logging
from bclib.obj_lib.lateral_obj import lateral_bc
import netCDF4 as nc
import numpy as np

class gib():
    def __init__(self,conf):
        self.config = conf
        self.gibilterra = lateral_bc(conf.input_data.file_nutrients)
    def nutrient_dataset_by_index(self,jn):
        if jn==0 : return self.gibilterra.phos
        if jn==1 : return self.gibilterra.ntra
        if jn==2 : return self.gibilterra.dox
        if jn==3 : return self.gibilterra.sica
        if jn==4 : return self.gibilterra.dic
        if jn==5 : return self.gibilterra.alk
        if jn==6 : return np.array()
    def generate(self,mask,bounmask_obj):
        if bounmask_obj.index is None:
            bounmask_obj.generate(mask)
        jpk, jpj, jpi = mask.shape
        vnudg = self.config.variables
        nudg = len(vnudg)
        
        for time in range(4):
            name_file = self.input_data.dir_out+"/GIB_yyyy"+self.gibilterra.season[time]+".nc"
            ncfile = nc.Dataset(name_file, 'w')
            for jn in range(nudg):
                GIB_matrix = self.nutrient_dataset_by_index(jn)
                aux = bounmask_obj.resto[jn,:,:,:]
                isNudg = (aux != 0) & (mask.mask)
                N_nudgingpoints = isNudg.sum()
                idx=np.zeros((N_nudgingpoints,),np.int32);
                data=np.zeros((N_nudgingpoints,),np.float64);
                count = 0
                for jk in range(jpk):
                    for jj in range(jpj):
                        for ji in range(jpi):
                            if (isNudg[jk,jj,ji]):
                                idx[count] = bounmask_obj.index[jk,jj,ji]
                                    data[count] = GIB_matrix[time,jk,jj,ji];
                                count = count+1
                                
    
    
                    ncfile.createDimension('gib_idxt_N1p',count)
                    idx_n1p = ncfile.createVariable('gib_idxt_N1p', 'i4', ('gib_idxt_N1p',))
                    idx_n1p[:] = idx[:]
                    data_n1p = ncfile.createVariable('gib_N1p', 'f',('gib_idxt_N1p',))
                    data_n1p[:] = data[:]
    
                if jn == 1:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = self.gibilterra.ntra[time,cord[0],cord[1],cord[2]];
                        count = count +1;
    
                    ncfile.createDimension('gib_idxt_N3n',count)
                    idx_n3n = ncfile.createVariable('gib_idxt_N3n', 'i', ('gib_idxt_N3n',))
                    idx_n3n[:] = idx[:]
                    data_n3n = ncfile.createVariable('gib_N3n', 'f',('gib_idxt_N3n',))
                    data_n3n[:] = data[:]
    
                if jn == 2:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = self.gibilterra.dox[time,cord[0],cord[1],cord[2]];
                        count = count +1;
    
                    ncfile.createDimension('gib_idxt_O2o',count)
                    idx_o2o = ncfile.createVariable('gib_idxt_O2o', 'i', ('gib_idxt_O2o',))
                    idx_o2o[:] = idx[:]
                    data_o2o = ncfile.createVariable('gib_O2o', 'f',('gib_idxt_O2o',))
                    data_o2o[:] = data[:]
    
                if jn == 3:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = self.gibilterra.sica[time,cord[0],cord[1],cord[2]];
                        count = count +1;
    
    
                    ncfile.createDimension('gib_idxt_N5s',count)
                    idx_n3n = ncfile.createVariable('gib_idxt_N5s', 'i', ('gib_idxt_N5s',))
                    idx_n3n[:] = idx[:]
                    data_n5s = ncfile.createVariable('gib_N5s', 'f',('gib_idxt_N5s',))
                    data_n5s[:] = data[:]
    
                if jn == 4:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = self.gibilterra.dic[time,cord[0],cord[1],cord[2]];
                        count = count +1;
    
                    ncfile.createDimension('gib_idxt_O3c',count)
                    idx_o3c = ncfile.createVariable('gib_idxt_O3c', 'i', ('gib_idxt_O3c',))
                    idx_o3c[:] = idx[:]
                    data_o3c = ncfile.createVariable('gib_O3c', 'f',('gib_idxt_O3c',))
                    data_o3c[:] = data[:]
    
                if jn == 5:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = self.gibilterra.alk[time,cord[0],cord[1],cord[2]];
                        count = count +1;
    
                    ncfile.createDimension('gib_idxt_O3h',count)
                    idx_o3h = ncfile.createVariable('gib_idxt_O3h', 'i', ('gib_idxt_O3h',))
                    idx_o3h[:] = idx[:]
                    data_o3h = ncfile.createVariable('gib_O3h', 'f',('gib_idxt_O3h',))
                    data_o3h[:] = data[:]
    
                if jn == 6:
                    for cord in isNudg[jn]:
                        idx[count] = index[cord[0],cord[1],cord[2]];
                        data[count] = 0.0025;
                        count = count +1;
    
                    ncfile.createDimension('gib_idxt_N6r',count)
                    idx_n6r = ncfile.createVariable('gib_idxt_N6r', 'i', ('gib_idxt_N6r',))
                    idx_n6r[:] = idx[:]
                    data_n6r = ncfile.createVariable('gib_N6r', 'f',('gib_idxt_N6r',))
                    data_n6r[:] = data[:]
    
    
                ncfile.close()