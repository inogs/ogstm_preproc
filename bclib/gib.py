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

    def generate(self, mask, bounmask_obj):
        
        logging.info("GIB files generation: start")
        
        if bounmask_obj.idx is None:
            bounmask_obj.generate(mask)
        
        # Common variables and definitions:
        jpk, jpj, jpi = mask.shape
        vnudg = self.config.variables
        nudg = len(vnudg)
        GIB_missing_values = np.ones((jpk, jpj, jpi), np.float64) * 1.e+20
        
        for time in range(4):
            
            # netCDF file preparation
            filename = self.config.dir_out + "GIB_yyyy" + self.gibilterra.season[time] + ".nc"
            ncfile = nc.Dataset(filename, 'w')
            ncfile.createDimension("lon", jpi)
            ncfile.createDimension("lat", jpj)
            ncfile.createDimension("dep", jpk)
            
            for jn in range(nudg):
                
                GIB_matrix = self.nutrient_dataset_by_index(jn)
                aux = bounmask_obj.resto[jn, ...]
                isNudg = (aux != 0) & (mask.mask)
                
                # Final dataset
                GIB_ready = np.where(isNudg, GIB_matrix[time, ...], GIB_missing_values)
                
                # Dump netCDF file
                vardataname = "gib_" + vnudg[jn][0]
                ncvar = ncfile.createVariable(vardataname, 'f', ("dep", "lat", "lon"))
                ncvar = GIB_ready
                
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
