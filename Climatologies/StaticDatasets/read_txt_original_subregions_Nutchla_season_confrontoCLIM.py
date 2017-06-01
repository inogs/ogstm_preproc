import numpy as np
import scipy.io.netcdf as NC
import pylab as pl
import glob
import os,sys
import time,datetime
THISDIR="/pico/home/userexternal/ssalon00/MS-REAN"
from maskload import *
import csv
#import matplotlib
#matplotlib.rc('text', usetex = True)

pl.close('all')

tk_m     = getDepthIndex(nav_lev,0.5)
mask_2D   = tmask[tk_m,:,:]
mask_2D[Lon < -5.3]   = False
topoNOGIB=mask_2D.copy()

IDrun='SO-05' 
FIGDIR=THISDIR
CLIMDIR='/fermi/home/userexternal/ssalon00/CLIMA_DB_Manca/OUTPUT_NEW/'
CLIMDIR2='/fermi/home/userexternal/ssalon00/CLIMA_DB_Leslie/'
MODDIR='/gpfs/scratch/userexternal/ssalon00/'+IDrun+'/wrkdir/POSTPROC/output/STAT_PROFILES/'

subreg_Manca_all=['DF2','DF3','DF4','DS2','DF1','DS1','DS3','DS4','DT1','DT2','DT3','DI1','DI3','DJ1','DJ2','DJ3','DJ6','DJ4','DJ7','DJ8','DJ5','DH1','DH2','DH3','DL1','DL2','DL3','DL4']
#subreg_Manca=['DF2','DT1','DS4','DS1','DJ5','DL1','DL3','DJ3','DI3']
#subreg_rLesl=['NWM','TYR','SWE','ALB','ION','LEV','LEV','ADS','ION'] #corrispondent
#subreg_Manca=['NWM','TYR','SWE','ALB','ION','LEV','LEV','ADS','ION']
#subreg_rLesl=['NWM','TYR','SWE','ALB','ION','LEV','LEV','ADS','ION'] #corrispondent
##subreg_Manca=['nwm','tyr','swe','alb','ion','lev','ads']
##subreg_rLesl=['NWM','TYR','SWE','ALB','ION','LEV','ADS'] #corrispondent
##subreg_in_fig=['NWM','TYR','SWE','ALB','ION','LEV','ADS']
## solo 2 alla volta per annual:
subreg_Manca=['lev','ads']
subreg_rLesl=['LEV','ADS'] #corrispondent
subreg_in_fig=['LEV','ADS']
#subreg_rLesl=['ALB','SWW','SWE','NWM','TYR','ADN','ADS','AEG','ION','LEV']
#subreg_mod=['alb','sww','swe','nwm','tyr','adn','ads','aeg','ion','lev']

var_ref_Manca=['CHL' ,'NO3' ,'PO4' ,'OXY' ,'SIL']
var_ref_rLesl=['CPHL','NTRA','PHOS','DOXY','SLCA']
var_mod_Manca=['P_i' ,'N3n' ,'N1p' ,'O2o' ,'N5s']

time_ref_rLesl=[    '',  '1-',  '2-',  '3-',  '4-']
time_ref_Manca=['annual','winter','spring','summer','autumn']
time_ref_in_fig=['annual','WINTER','SPRING','SUMMER','AUTUMN']

xmin=[0, 0,0,  180,0]
xmax=[1,10,0.6,280,9]

TEXTxlabel3=['mgChl/m'+u'\u00B3','mmolN/m'+u'\u00B3','mmolP/m'+u'\u00B3','mmolO'+u'\u2082'+'/m'+u'\u00B3','mmolSi/m'+u'\u00B3']

nt=int(sys.argv[1])
if(nt==0): # annual statistics
   m0=0 
if(nt==1): # winter statistics JFM
   m0=0 # Jan
if(nt==2): # spring statistics AMJ
   m0=3 # Apr
if(nt==3): # summer statistics JAS
   m0=6 # Jul
if(nt==4): # autumn statistics OND
   m0=9 # Oct

#nt=0 # only annual statistics
print(nt,time_ref_Manca[nt])

in_year=1961
fn_year=1999
print('postproc from '+str(in_year)+' to '+str(fn_year))
nyears=fn_year-in_year+1

# define dates available from model output
alldates=[]
dd=[]
datelist=glob.glob(MODDIR + "ave.19*profiles.nc")
datelist.sort()
for date in datelist:
    data=os.path.basename(date)
    alldates.append(data)
    dd.append(data[4:12])

ntot=len(dd)


for ns in range(len(subreg_Manca)):
#for ns in range(0,len(subreg_Manca)/2):
#for ns in range(len(subreg_Manca)/2,len(subreg_Manca)):

    fig1=pl.figure(num=None,figsize=(9,6),dpi=300,facecolor='w',edgecolor='k')
    ax0=fig1.add_subplot(2,len(var_ref_Manca)+1,1+len(var_ref_Manca)+1)
    M=NC.netcdf_file(submaskfile,"r")
    s=M.variables[subreg_Manca[ns]].data[0,:,:]
    M.close()
    s0=s.copy()
    s0[s0==0]=np.NaN
    pl.imshow(s0[:,55:])
    pl.contour(topoNOGIB[:,55:],[0],linewidths=0.5,colors='k')
    pl.gca().invert_yaxis()
    ax0.set_xticklabels([])
    ax0.set_yticklabels([])

    for nv in range(len(var_ref_Manca)):

# READING CLIMATOLOGICAL REF DATA Manca et al. (2004)
        filename=var_ref_Manca[nv] + '/' + time_ref_Manca[nt] + '.' + subreg_Manca[ns] + '.txt'
        print(filename)
        n=0
        with open(CLIMDIR + filename,'rb') as f:
            reader=csv.reader(f)
            for row in reader:
                n=n+1

        if (n>0): # if data exist read them  
             print('n.levels=',n-1)

             climm=np.zeros((n,4)) #z,mean,std,number
             head=[]

             crs = open(CLIMDIR + filename,"r")
             i=0
             for columns in ( raw.strip().split() for raw in crs ):
                if(i==0):
                   head=columns
                else:
                   climm[i,0]=columns[0] #zlev
                   climm[i,1]=columns[1] #mean
                   climm[i,2]=columns[2] #std
                   climm[i,3]=columns[3] #num.points
# convert dissolved oxygen:
                   if (nv==3):
                       climm[i,1:3]=climm[i,1:3]*44.66
#http://guillaumemaze.googlecode.com/svn/trunk/matlab/codes/geophysic/convert_oxygen.m
                i=i+1

             clim=np.zeros((n-1,4)) #z,mean,std,number
             clim[:,:]=climm[1:,:]
        else:   # if data not exist create a fake clim matrix
             clim=np.zeros((1,4))*np.NaN
#
# READING CLIMATOLOGICAL REF DATA Leslie (2013)
        filename2='plot-average-' + time_ref_rLesl[nt] + subreg_rLesl[ns] + '-' + var_ref_rLesl[nv] + '.csv'
    #    try:
    #        with open('filename'): pass
    #    except IOError:
    #        nv=nv+1
    #        filename='plot-median-' + time_ref[nt] + subreg_ref[ns] + '-' + var_ref[nv] + '.csv'
    #
        print(filename2)
        n=0
        with open(CLIMDIR2 + filename2,'rb') as f:
            reader=csv.reader(f)
            for row in reader:
                n=n+1

        print('n.lines=',n)

        clim2=np.zeros((n,4)) #z,mean,std,number

        quotes=open(CLIMDIR2 + filename2, "rU" )
        #titles= quotes.next().strip().split( ',' )
        i=0
        for q in quotes:
            values= q.strip().split( ',' )
            #data= dict( zip(titles,values) )
            #print values
            clim2[i,0]=float(values[1]) #zlev
            clim2[i,1]=float(values[0]) #mean
            clim2[i,2]=float(values[2]) #std
            clim2[i,3]=float(values[3]) #num.points
# convert dissolved oxygen:
            if (nv==3):
                clim2[i,1:4]=clim2[i,1:4]*(44.66/1.42903)
#http://guillaumemaze.googlecode.com/svn/trunk/matlab/codes/geophysic/convert_oxygen.m
            i=i+1
# READING MODEL DATA (Mean + STD)
        # nc=1 area with depth>200m (off-shore)
        nc=1 #always this coherently with climatology reference

######## Fig.1 ALL DATA
######## down to 200m
        ax1=fig1.add_subplot(2,len(var_ref_Manca),nv+1)
        for ny in range(nyears): #  PLOT one red line for each month ... of the four given seasons
            m1=m0+ny*12
            m2=m1+3
            if(nt==0):
               m2=m1+12
            nm=0
            #print('year:'+str(ny))
            for infile in datelist[m1:m2]:  # read 3 month of a given season in the given year
                F=NC.netcdf_file(infile,"r")
                sublist_F=F.sub___list.replace(' ','').rsplit(',')
                indsub=sublist_F.index(subreg_Manca[ns]) # trovo l'indice del subbasin dalla lista delle subbasin salvae
                t=F.variables[var_mod_Manca[nv]].data[indsub,nc,:,:]
                F.close()
                t0=t.copy()
            #    pl.plot(t0[:,0],-nav_lev[:],'-r')

        for ny in range(nyears): # PLOT MEAN OF THE SEASON FOR A GIVEN YEAR
            m1=m0+ny*12
            m2=m1+3
            tm=np.zeros((len(nav_lev),5,3))
            if(nt==0):
               m2=m1+12
               tm=np.zeros((len(nav_lev),5,12))
            nm=0
            #print('year:'+str(ny))
            for infile in datelist[m1:m2]:
                F=NC.netcdf_file(infile,"r")
                sublist_F=F.sub___list.replace(' ','').rsplit(',')
                indsub=sublist_F.index(subreg_Manca[ns]) # trovo l'indice del subbasin dalla lista delle subbasin salvae
                #print(infile)
                tm[:,:,nm]=F.variables[var_mod_Manca[nv]].data[indsub,nc,:,:]
                F.close()
                nm=nm+1
            color=str(1.-float(ny)/nyears) #to create greyscale for all the years plotted
            #print(ny,color)
            ii=tm[:,0,0]>0.
            pl.plot(tm[:,0,:].mean(axis=1),-nav_lev[:],color,label=str(in_year+ny))
            #if(nv==0):
            #    pl.plot(tm[:,0].mean(axis=1),-nav_lev[:],color,label=str(in_year+ny))
            #else:
            #    pl.plot(tm[ii,0].mean(axis=1),-nav_lev[ii],color,label=str(in_year+ny))

        pl.plot(clim[:,1],-clim[:,0],'.b') # PLOT CLIM MANCA
        pl.plot(clim[:,1]-clim[:,2],-clim[:,0],'-.b')
        pl.plot(clim[:,1]+clim[:,2],-clim[:,0],'-.b')
#        pl.plot(clim2[:,1],clim2[:,0],'.b') # PLOT CLIM LESLIE
#        pl.plot(clim2[:,1]-clim2[:,2],clim2[:,0],'-.b')
#        pl.plot(clim2[:,1]+clim2[:,2],clim2[:,0],'-.b')
#        ax1.set_title(subreg_Manca[ns] + ' ' + var_ref[nv] + ' (median)')
        ax1.set_title(var_ref_Manca[nv])
        ax1.grid(color='k',linestyle='--')
        ax1.set_xlabel(TEXTxlabel3[nv])
#        ax1.set_ylabel('m')
        if(ns==3):
           xmax[0]=1.5
        else:
           xmax[0]=1.0
        pl.xlim(xmin[nv],xmax[nv])
        pl.ylim(-200,0)
        pl.rc('xtick', labelsize=8)
        if(nv==0):
            pl.rc('ytick', labelsize=8)
            ax1.set_xlabel(TEXTxlabel3[nv])
        else:
            ax1.set_yticklabels([])           

        if (nv==4):
            ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,labelspacing=-0.25, handletextpad=0,borderpad=0.1)
            leg = pl.gca().get_legend()
            ltext  = leg.get_texts()
            pl.setp(ltext,fontsize=8)

######## from 200m to 2000m
        if (nv>0):
            ax2=fig1.add_subplot(2,len(var_ref_Manca),nv+1+len(var_ref_Manca))
            for ny in range(nyears): #  PLOT one red line for each month ... of the four given seasons
                m1=m0+ny*12
                m2=m1+3
                if(nt==0):
                   m2=m1+12
                nm=0
            #print('year:'+str(ny))
                for infile in datelist[m1:m2]:  # read 3 month of a given season in the given year
                    F=NC.netcdf_file(infile,"r")
                    sublist_F=F.sub___list.replace(' ','').rsplit(',')
                    indsub=sublist_F.index(subreg_Manca[ns]) # trovo l'indice del subbasin dalla lista delle subbasin salvae
                    t=F.variables[var_mod_Manca[nv]].data[indsub,nc,:,:]
                    F.close()
                    t0=t.copy()
            #        pl.plot(t0[:,0],-nav_lev[:],'-r')

            for ny in range(nyears): # PLOT MEAN OF THE SEASON FOR A GIVEN YEAR
                m1=m0+ny*12
                m2=m1+3
                tm=np.zeros((len(nav_lev),5,3))
                if(nt==0):
                   m2=m1+12
                   tm=np.zeros((len(nav_lev),5,12))
                nm=0
                for infile in datelist[m1:m2]:
                    F=NC.netcdf_file(infile,"r")
                    sublist_F=F.sub___list.replace(' ','').rsplit(',')
                    indsub=sublist_F.index(subreg_Manca[ns]) # trovo l'indice del subbasin dalla lista delle subbasin salvae
                    tm[:,:,nm]=F.variables[var_mod_Manca[nv]].data[indsub,nc,:,:]
                    F.close()
                    nm=nm+1
                color=str(1.-float(ny)/nyears) #to create greyscale for all the years plotted
                #print(ny,color)
                ii=tm[:,0,0]>0.
                lastlev=nav_lev[ii][-1]
                pl.plot(tm[:,0,:].mean(axis=1),-nav_lev[:],color,label=str(in_year+ny))
                #if(nv==0):
                #    pl.plot(tm[:lx,0].mean(axis=1),-nav_lev[:lx],color,label=str(in_year+ny))
                #else:
                #    pl.plot(tm[ii,0].mean(axis=1),-nav_lev[ii],color,label=str(in_year+ny))


            pl.plot(clim[:,1],-clim[:,0],'.b') # PLOT CLIM MANCA
            pl.plot(clim[:,1]-clim[:,2],-clim[:,0],'-.b')
            pl.plot(clim[:,1]+clim[:,2],-clim[:,0],'-.b')
#            pl.plot(clim2[:,1],clim2[:,0],'.b') # PLOT CLIM LESLIE
#            pl.plot(clim2[:,1]-clim2[:,2],clim2[:,0],'-.b')
#            pl.plot(clim2[:,1]+clim2[:,2],clim2[:,0],'-.b')
#            ax1.set_title(subreg_Manca[ns] + ' ' + var_ref[nv] + ' (median)')
#            ax2.set_title(var_ref_Manca[nv])
            ax2.grid(color='k',linestyle='--')
            ax2.set_xlabel(TEXTxlabel3[nv])
            pl.xlim(xmin[nv],xmax[nv])
            pl.ylim(-lastlev,-200)           
#            ax1.set_ylabel('m')
            pl.rc('xtick', labelsize=8)
            if (nv==1):
                pl.rc('ytick', labelsize=8)
            else:
                ax2.set_yticklabels([]) 

        #pl.show(mainloop=False)

    fig1.text(0.12,0.40,time_ref_in_fig[nt],fontsize=15)
    fig1.text(0.12,0.35,IDrun,fontsize=15)
    fig1.text(0.12,0.15,subreg_in_fig[ns],fontsize=15)
    fig1.savefig(FIGDIR + 'prof_mean_confronto_' + time_ref_Manca[nt] + '_' + subreg_Manca[ns] + '.png',dpi=fig1.dpi)

print('EOB')
