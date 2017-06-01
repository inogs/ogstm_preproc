import scipy.io.netcdf as NC
from config import LayerList, REQUESTORS_LIST, basV2

climfile="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Nutrients/clim.phosphate.nc"
ncIN=NC.netcdf_file(climfile,'r')
PHOS = ncIN.variables['CLIMATOLOGY'].data.copy()
ncIN.close()

climfile="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Nutrients/clim.nitrate.nc"
ncIN=NC.netcdf_file(climfile,'r')
NIT = ncIN.variables['CLIMATOLOGY'].data.copy()
ncIN.close()

climfile="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Nutrients/clim.oxygen.nc"
ncIN=NC.netcdf_file(climfile,'r')
DOX = ncIN.variables['CLIMATOLOGY'].data.copy()
ncIN.close()

climfile="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Carbonatics/clim.ALK.nc"
ncIN=NC.netcdf_file(climfile,'r')
ALK = ncIN.variables['CLIMATOLOGY'].data.copy()
ncIN.close()

climfile="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/Carbonatics/clim.DIC.nc"
ncIN=NC.netcdf_file(climfile,'r')
DIC = ncIN.variables['CLIMATOLOGY'].data.copy()
ncIN.close()


StatList=['Mean','Std','p01','p25','p50','p75','p99','MeanDepth','MeanMonth','nValues']

def statistics(varname, iSeas,iSub):
    '''
    Returns:
    * STAT * (depths=14,stats=10) numpy array
    '''
    if varname == "N1p" : return PHOS[iSeas,iSub,:,:]
    if varname == "O2o" : return  DOX[iSeas,iSub,:,:]
    if varname == "N3n" : return  NIT[iSeas,iSub,:,:]
    if varname == "Ac"  : return  ALK[iSeas,iSub,:,:]
    if varname == "DIC" : return  DIC[iSeas,iSub,:,:]



if __name__=="__main__":
    for iSeas, SeasReq in enumerate(REQUESTORS_LIST):
        for iSub, sub in enumerate(basV2.P):
            for ilayer, layer in enumerate(LayerList):
                STAT = statistics('N1p',iSeas, iSub)