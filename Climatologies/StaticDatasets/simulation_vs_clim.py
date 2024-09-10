from config import LayerList, REQUESTORS_LIST, basV2
import clim_reader
import figure_generator
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.timeseries.plot import Hovmoeller_matrix
import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
import pylab as pl
from bitsea.commons.utils import getcolor
IDrun='eas_11'
OUTDIR="IMG/"
MODDIR="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
TI = TimeInterval("20140101","20180101","%Y%m%d")
maskfile ="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc"
maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
Mask8 = Mask(maskfile8)
TheMask= Mask(maskfile)
jpk,jpj,jpi = TheMask.shape
z = -TheMask.zlevels

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")

REQUESTORS_FOR_PLOTLINE=TL.getMonthlist()
nPlotlines=len(REQUESTORS_FOR_PLOTLINE)
def getcolor2(ntimes,itime):
    '''
    Arguments:
    * nTimes * integer, the number of times
    * itime  * integer, time index

    if n. years <  8: use different colors
    if n. years >= 8: use grey scale color

    Returns :
    * color * string
    '''
    yearcol = ['k','b','r','g','c','m','y']
    if ntimes < 8:
        color=yearcol[itime]
    else:
        color=str(1.-float(itime+1)/ntimes)
    return color


# SECTION FIGURE 3

def TimelistSelection(ii_seas, ii):
    '''
    Arguments
    * ii_seas * list of indexes of that season
    * ii      * list of indexes (generic)
    '''
    datetimelist= [TL.Timelist[k] for k in ii_seas if k in ii]
    filelist    = [TL.filelist[k] for k in ii_seas if k in ii]
    datetimelist= [TL.Timelist[k] for k in ii]
    filelist    = [TL.filelist[k] for k in ii]
    return datetimelist, filelist

VARLIST=['P_l','N1p','N3n','O2o','N5s']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

for iSeas, SeasReq in enumerate(REQUESTORS_LIST):
    if iSeas<4 : continue
    ii_seas,w = TL.select(SeasReq)
    for iSub, sub in enumerate(basV2.P):
        submask = SubMask(sub,maskobject=Mask8)
        F = figure_generator.figure_generator(submask)
        fig, axes = F.gen_structure_1(IDrun,SeasReq.string,sub.name)
        outfile = OUTDIR + "prof_mean_confronto_1_" + SeasReq.string + "." + sub.name + ".png"
        
        for iTime, req in enumerate(REQUESTORS_FOR_PLOTLINE):
            ii,w = TL.select(req)
            datetimelist, filelist = TimelistSelection(ii_seas, ii)

            MEAN = np.zeros((jpk,), dtype=var_dype )
            for ivar, var in enumerate(VARLIST):
                Mean_profiles,_,_ = Hovmoeller_matrix(datetimelist,filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
                mean_profile = Mean_profiles.mean(axis=1)
                mean_profile[mean_profile==0]=np.nan
                MEAN[var] = mean_profile
            
            label = req.string
            color = getcolor(nPlotlines, iTime)
            figure_generator.profile_plotter(z,MEAN['P_l'],color, axes[0], None,   label)
            figure_generator.profile_plotter(z,MEAN['N3n'],color, axes[1], axes[6],label)
            figure_generator.profile_plotter(z,MEAN['N1p'],color, axes[2], axes[7],label)
            figure_generator.profile_plotter(z,MEAN['O2o'],color, axes[3], axes[8],label)
            figure_generator.profile_plotter(z,MEAN['N5s'],color, axes[4], axes[9],label)            
        
        N3n_STAT = clim_reader.statistics('N3n', iSeas, iSub)
        N1p_STAT = clim_reader.statistics('N1p', iSeas, iSub)
        O2o_STAT = clim_reader.statistics('O2o', iSeas, iSub)
        N3n_clim_mean = N3n_STAT[:,0]
        N1p_clim_mean = N1p_STAT[:,0]
        O2o_clim_mean = O2o_STAT[:,0]

        N3n_clim_std  = N3n_STAT[:,1]
        N1p_clim_std  = N1p_STAT[:,1]
        O2o_clim_std  = O2o_STAT[:,1]
        
        figure_generator.clim_profile_plotter(z_clim,N3n_clim_mean,N3n_clim_std, axes[1], axes[6])
        figure_generator.clim_profile_plotter(z_clim,N1p_clim_mean,N1p_clim_std, axes[2], axes[7])
        figure_generator.clim_profile_plotter(z_clim,O2o_clim_mean,O2o_clim_std, axes[3], axes[8])

        figure_generator.add_legend(axes)
        fig.savefig(outfile)
        print outfile
        pl.close(fig)



# SECTION FIGURE 2

VARLIST=['P_l','P_c','ppn','B1c']
# per ora non mettiamo la PAR
#VARLIST=['P_l','PAR','P_c','ppn','B1c']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)

for iSeas, SeasReq in enumerate(REQUESTORS_LIST):
    if iSeas<4 : continue
    ii_seas,w = TL.select(SeasReq)
    for iSub, sub in enumerate(basV2.P):
        submask = SubMask(sub,maskobject=Mask8)
        F = figure_generator.figure_generator(submask)
        fig, axes = F.gen_structure_2(IDrun,SeasReq.string,sub.name)
        outfile = OUTDIR + "prof_mean_confronto_2_" + SeasReq.string + "." + sub.name + ".png"
        
        for iTime, req in enumerate(REQUESTORS_FOR_PLOTLINE):
            ii,w = TL.select(req)
            datetimelist, filelist = TimelistSelection(ii_seas, ii)

            MEAN = np.zeros((jpk,), dtype=var_dype )
            for ivar, var in enumerate(VARLIST):
                Mean_profiles,_,_ = Hovmoeller_matrix(datetimelist,filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
                mean_profile = Mean_profiles.mean(axis=1)
                mean_profile[mean_profile==0]=np.nan
                MEAN[var] = mean_profile

            label = req.string
            color = getcolor(nPlotlines, iTime)
            figure_generator.profile_plotter(z,MEAN['P_l'],color, axes[0], None,   label)
            #figure_generator.profile_plotter(z,MEAN['PAR'],color, axes[1], axes[6],label)
            figure_generator.profile_plotter(z,MEAN['P_c'],color, axes[2], axes[7],label)
            figure_generator.profile_plotter(z,MEAN['ppn'],color, axes[3], axes[8],label)
            figure_generator.profile_plotter(z,MEAN['B1c'],color, axes[4], axes[9],label)
            
        figure_generator.add_legend(axes)
        fig.savefig(outfile)
        print outfile
        pl.close(fig)



# SECTION FIGURE 3

VARLIST=['pCO2','DIC','Ac','pH','R2c']
var_dype = [(var,np.float32) for var in VARLIST]
nVar = len(VARLIST)


for iSeas, SeasReq in enumerate(REQUESTORS_LIST):
    if iSeas<4 : continue
    ii_seas,w = TL.select(SeasReq)
    for iSub, sub in enumerate(basV2.P):
        submask = SubMask(sub,maskobject=Mask8)
        F = figure_generator.figure_generator(submask)
        fig, axes = F.gen_structure_3(IDrun,SeasReq.string,sub.name)
        outfile = OUTDIR + "prof_mean_confronto_3_" + SeasReq.string + "." + sub.name + ".png"
        
        for iTime, req in enumerate(REQUESTORS_FOR_PLOTLINE):
            ii,w = TL.select(req)
            datetimelist, filelist = TimelistSelection(ii_seas, ii)

            MEAN = np.zeros((jpk,), dtype=var_dype )
            for ivar, var in enumerate(VARLIST):
                Mean_profiles,_,_ = Hovmoeller_matrix(datetimelist,filelist, var, iSub, coast=1, stat=0, depths=np.arange(jpk)) #72 nFiles
                mean_profile = Mean_profiles.mean(axis=1)
                mean_profile[mean_profile==0]=np.nan
                MEAN[var] = mean_profile
            

            label = req.string
            color = getcolor(nPlotlines, iTime)
            figure_generator.profile_plotter(z,MEAN['pCO2'],color, axes[0], None,   label)
            figure_generator.profile_plotter(z,MEAN['DIC'],color, axes[1], axes[6],label)
            figure_generator.profile_plotter(z,MEAN['Ac'],color, axes[2], axes[7],label)
            figure_generator.profile_plotter(z,MEAN['pH'],color, axes[3], axes[8],label)
            figure_generator.profile_plotter(z,MEAN['R2c'],color, axes[4], axes[9],label)
            
        
        DIC_STAT = clim_reader.statistics('DIC', iSeas, iSub)
        Ac__STAT = clim_reader.statistics('Ac', iSeas, iSub)

        Dic_clim_mean = DIC_STAT[:,0]
        Ac__clim_mean = Ac__STAT[:,0]


        Dic_clim_std  = DIC_STAT[:,1]
        Ac__clim_std  = Ac__STAT[:,1]

        
        figure_generator.clim_profile_plotter(z_clim,Dic_clim_mean,Dic_clim_std, axes[1], axes[6])
        figure_generator.clim_profile_plotter(z_clim,Ac__clim_mean,Ac__clim_std, axes[2], axes[7])

        figure_generator.add_legend(axes)
        fig.savefig(outfile)
        print outfile
        pl.close(fig)





