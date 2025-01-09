import logging
from bclib.lateral import lateral_bc
import netCDF4 as nc
import numpy as np
import glob, os

class gib():
    def __init__(self,conf,mask):
        self.config = conf
        logging.info("Lateral conditions : start")
        self.gibilterra = lateral_bc(conf,mask)
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
        if jn==7 : return self.gibilterra.Hg0
        if jn==8 : return self.gibilterra.Hg2
        if jn==9 : return self.gibilterra.MMHg
        if jn==10: return self.gibilterra.DMHg

    def generate(self, mask, bounmask_obj, all_variables=True):
        
        logging.info("GIB files generation: start")
        
        if bounmask_obj.idx is None:
            bounmask_obj.generate(mask)
        
        # Common variables and definitions:
        jpk, jpj, jpi = mask.shape
        vnudg = self.config.variables
        nudg = len(vnudg)
        nudg_classic_variables=[ vnudg[k][0] for k in range(nudg)]
        RST_LIST=glob.glob(self.config.RST_FILES)
        OTHER_VARIABLES=[]
        for filename in RST_LIST:
            basename = os.path.basename(filename)
            var  = basename.rsplit(".")[2]
            if var not in nudg_classic_variables:
                OTHER_VARIABLES.append(var)
        #GIB_missing_values = np.ones((jpk, jpj, jpi), np.float64) * 1.e+20
        jpt = len(self.config.gib_season)
        for time in range(jpt):
            
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
                print (jn, nudg, vnudg[jn][0])
                
                GIB_matrix[time, ~isNudg] = -1.
                GIB_matrix[time, ~mask.mask] = 1.e+20
                
                # Dump netCDF file
                ncvar = ncfile.createVariable(vnudg[jn][0], 'f', ("dep", "lat", "lon"))
                setattr(ncvar, 'missing_value', np.float32(1.e+20))
                ncvar[...] = GIB_matrix[time, ...]

            aux = bounmask_obj.resto[0, ...] # the same of N1p
            isNudg = (aux != 0) & (mask.mask)
            for filename in RST_LIST:
                basename = os.path.basename(filename)
                var  = basename.rsplit(".")[2]
                if not all_variables : continue
                if var in nudg_classic_variables:  continue
                GIB_matrix = self.read(filename, "TRN" + var)
                GIB_matrix[0, ~isNudg] = -1.
                GIB_matrix[0, ~mask.mask] = 1.e+20
                ncvar = ncfile.createVariable(var, 'f', ("dep", "lat", "lon"))
                ncvar[:] = GIB_matrix[0, ...]
                setattr(ncvar, 'missing_value', np.float32(1.e+20))
                
            ncfile.close()
            
        logging.info("GIB files generation: done")
        
if __name__ == "__main__":
    from bitsea.commons.mask import Mask
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
