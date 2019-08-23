from commons import netcdf4
from commons.Timelist import TimeList,TimeInterval
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import numpy as np
from commons import timerequestors 

TheMask=Mask("/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_02/wrkdir/MODEL/meshmask.nc")
Mask_0 = TheMask.cut_at_level(0)
area = TheMask.area
jpk, jpj, jpi= TheMask.shape

INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_02/wrkdir/MODEL/BC"
TI = TimeInterval("2017","2018","%Y")
TL =TimeList.fromfilenames(TI, INPUTDIR, "TIN_2*", prefix="TIN_")
nFrames = TL.nTimes

bool_dtype=[(sub.name,np.bool)    for sub in OGS.P]
sub__dtype=[(sub.name,np.float32) for sub in OGS.P]

SUBM = np.zeros((jpj, jpi), dtype=bool_dtype)
for isub, sub in enumerate(OGS.Pred):
    S=SubMask(sub, maskobject=Mask_0)
    SUBM[sub.name]=S.mask
    SUBM['med']=(SUBM['med'] | SUBM[sub.name])

nSub = len(OGS.P.basin_list)

def conversion(var):
    '''
    from mmol or mg to KTONS(N1p,N3n,N5s,O3c,O3h)  or Gmol (O2o)
    '''
    Conversion={}
    w= 1.0e-12
    n = 14
    p = 31
    s = 28
    Conversion['N1p'] =w*p
    Conversion['N3n'] =w*n
    Conversion['N5s'] =w*s
    Conversion['O3h'] =w
    Conversion['O3c'] =w
    Conversion['O2o'] =w
    return Conversion[var]

VARLIST=["N3n","N1p","O2o","N5s","O3h","O3c"]
BALANCE_KT_MONTH={}
for var in VARLIST:
    print var
    MMOL_MONTH = np.zeros((nFrames,),dtype=sub__dtype)
    KTON_MONTH = np.zeros((nFrames,),dtype=sub__dtype)
    for iFrame, filename in enumerate(TL.filelist):
        req = timerequestors.Monthly_req(2017,iFrame+1)
        B=netcdf4.readfile(filename, "riv_" + var)
        good = B > -1
        for isub, sub in enumerate(OGS.P):
            localmask=SUBM[sub.name]
            V = B*area
            river_on_sub = V[localmask & good].sum() #mmol/s or mg/s
            time = req.time_interval.length() 
            MMOL_MONTH[sub.name][iFrame] = river_on_sub *time #mmol/month
            KTON_MONTH[sub.name][iFrame] = river_on_sub *time * conversion(var)
    BALANCE_KT_MONTH[var]=KTON_MONTH

print BALANCE_KT_MONTH['N3n']['adr1']
print BALANCE_KT_MONTH['N3n']['adr1'].sum()

nvars = len(VARLIST)
TABLE= np.zeros((nSub,nvars),np.float32)
for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(OGS.P):
        TABLE[isub,ivar]= BALANCE_KT_MONTH[var][sub.name].sum()
from commons.utils import writetable
rows_names_list= [sub.name for sub in OGS.P]
writetable('KTon_2017_subbasins.txt', TABLE, rows_names_list, VARLIST)

    



import config as conf
from bclib.river import river
R = river(conf)
R.gen_map_indexes(TheMask)

LINES=[]
for sub in OGS.Pred:
    LINES.append("********** " +  sub.name +  "*************\n" )
    for iriver in range(R.nrivers):
        ji=R.georef['indLon'][iriver]-1
        jj=R.georef['indLat'][iriver]-1
        if SUBM[sub.name][jj,ji]: 
            line="  %d       %s " %(iriver+1, R.xls_data['monthly'][iriver+1,5])
            LINES.append(line+"\n")
fid=open('rivernames_on_subbasin.txt','w')
fid.writelines(LINES)
fid.close()

#N3n on xls
#PO 105.8
#Piave Tagliamento Isonzo Livenza Brenta Lika Reno Krka 60.37
#Neretva 4.29


