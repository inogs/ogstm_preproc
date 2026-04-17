import config as conf
from bclib.river import river
from bclib.bounmask import bounmask
from bitsea.commons.mask import Mask
import numpy as np

BOUN = bounmask(conf)
index_inv = BOUN.load('index_inv')
TheMask = Mask.from_file(conf.file_mask)


filename=conf.dir_out + "TIN_20040515-00:00:00.nc"
R=river(conf)
idxt = R.load_from_file(filename, 'riv_idxt')
N1p  = R.load_from_file(filename, 'riv_N1p' )



w= 1.0e+12;
t = 1./(365 * 86400)
n = 1./14;
p = 1./31;
s = 1./28;
cn = w*t*n
cp = w*t*p
cs = w*t*s

diff = np.zeros((R.nrivers,),np.float32)
for ip in range(R.nrivers):
    jk,jj, ji = index_inv[idxt[ip]-1,:]
    print jk, ji, jj # similar to xsl file
    diff[ip] = R.xls_data['DIP_KTperYR_NOBLS']['2004'][ip] * R.monthly_mod[ip,4]/100*12*cp/TheMask.area[jj-1,ji-1] - N1p[ip]
print diff.max()

V = R.load_from_file(filename, 'riv_N3n')
sheet = "DIN_KTperYR_NOBLS"
conv=cn

for ip in range(R.nrivers):
    jk,jj, ji = index_inv[idxt[ip]-1,:]
    diff[ip] = R.xls_data[sheet]['2004'][ip] * R.monthly_mod[ip,4]/100*12*conv/TheMask.area[jj-1,ji-1] - V[ip]
print diff.max()


from bitsea.layer_integral.mapplot import mapplot
Map2d = TheMask.mask_at_level(0)
fig, ax = mapplot({'data':Map2d, 'clim':[0,1]})
xp = R.forced_coord[:,0]-1
yp = R.forced_coord[:,1]-1
ax.plot(xp, yp,'w.')
fig.set_dpi(500)
#fig.savefig('pippo.png')
