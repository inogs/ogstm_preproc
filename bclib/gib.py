import logging
from bclib.lateral import lateral_bc
import netCDF4 as nc
import numpy as np

class gib():
    def __init__(self,conf,mask):
        self.config = conf
        logging.info("Lateral conditions : start")
        self.gibilterra = lateral_bc(conf.file_nutrients,mask)
        logging.info("Lateral conditions done")

    def read(self,filename, var):
        ''' Useful in check and debug'''
        ncfile = nc.Dataset(filename, 'r')
        M = np.array(ncfile.variables[var])
        ncfile.close()
        return M
    def nutrient_dataset_by_index(self,jn):
        if jn==0 : return self.gibilterra.phos
        if jn==1 : return self.gibilterra.ntra
        if jn==2 : return self.gibilterra.dox
        if jn==3 : return self.gibilterra.sica
        if jn==4 : return self.gibilterra.dic
        if jn==5 : return self.gibilterra.alk
        if jn==6 :
            size_nutrients = self.gibilterra.phos.shape
            return np.ones(size_nutrients,np.float64)*0.0025  # G. Cossarini estimate

    def generate(self,mask,bounmask_obj):
        logging.info("GIB files generation: start")
        if bounmask_obj.idx is None:
            bounmask_obj.generate(mask)
        jpk, jpj, jpi = mask.shape
        vnudg = self.config.variables
        nudg = len(vnudg)
        
        for time in range(4):
            filename = self.config.dir_out + "GIB_yyyy" + self.gibilterra.season[time]+".nc"
            print filename
            ncfile = nc.Dataset(filename, 'w')
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
                                idx[count] = bounmask_obj.idx[jk,jj,ji]
                                data[count] = GIB_matrix[time,jk,jj,ji]
                                count = count+1


                dimension_name = "gib_idxt_" + vnudg[jn][0]
                vardataname    = "gib_"      + vnudg[jn][0]
                ncfile.createDimension(dimension_name,count)
                ncvar = ncfile.createVariable(dimension_name, 'i4', (dimension_name,))
                ncvar[:] = idx[:]
                ncvar = ncfile.createVariable(vardataname, 'f',(dimension_name,))
                ncvar[:] = data[:]
            ncfile.close()
        logging.info("GIB files generation: done")
if __name__ == "__main__":
    from commons.mask import Mask
    import config as conf
    conf.file_nutrients = "../"+ conf.file_nutrients
    conf.file_mask="../masks/meshmask_872.nc"
    conf.dir_out = "../out"
    TheMask = Mask(conf.file_mask)
    GIB = gib(conf,TheMask)
    from bounmask import bounmask
    BOUN=bounmask(conf)
    GIB.generate(TheMask, BOUN)
    old_values = GIB.read("GIB_yyyy0215-12:00:00.nc","gib_N3n")
    new_values = GIB.read("../out/GIB_yyyy0215-12:00:00.nc","gib_N3n")

    import pylab as pl
    fig,ax=pl.subplots()
    d=new_values-old_values
    ax.plot(d)
    fig.show()
