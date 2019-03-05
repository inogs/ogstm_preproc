from commons.Timelist import TimeInterval, TimeList
import numpy as np
import pylab as pl
from commons.utils import writetable
#TI=TimeInterval("20161001","20161031","%Y%m%d")
INPUTDIR="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_07/wrkdir/MODEL/FORCINGS/"
TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", filtervar="U", prefix="U",hour=0)

outfile="/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/deltaT.txt"
VALUES=np.loadtxt(outfile)
nFrames,_=VALUES.shape
DELTAT=VALUES[:,0]
NUMBER=VALUES[:,1]
K = VALUES[:,2]
J = VALUES[:,3]
I = VALUES[:,4]

fig,ax=pl.subplots()

ax2 = ax.twinx()
ax2.bar(TL.Timelist,NUMBER, width=3, color='0.8', label='# deltat <450 s')

ax.plot(TL.Timelist, DELTAT,'r',label="Delta t")
ax.grid()
ax.set_ylabel('delta t (s)')
ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
ax.patch.set_visible(False) # hide the 'canvas'
ax.legend()
ax2.set_ylabel('# deltat < 450')
fig.show()
fig.set_size_inches(16.0, 4)
fig.set_dpi(150)
fig.savefig('deltat.png')

rows_names_list=[time.strftime("%Y%m%d") for time in TL.Timelist]
column_names_list=["deltat (s)", "# (deltat < 450)", "Kmax", "Jmax", "Imax"]
writetable("DELTAT", VALUES, rows_names_list, column_names_list, fmt="%10.3f  %d %d %d %d")

for i in range(nFrames):
    if NUMBER[i] > 0:
        if ((J[i]>70) & (J[i]<170) ):
            if ((I[i]> 200) & (I[i]<382)):
                print "%s %10.3f %d " % (TL.Timelist[i].strftime("%Y%m%d"), DELTAT[i], NUMBER[i] )
