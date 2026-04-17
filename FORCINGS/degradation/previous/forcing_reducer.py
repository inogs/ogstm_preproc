import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Degradates forcings INGV files in another resolution defined for ogstm model.
    New resolution can be 2,3,4,... coarser than the original one.
    In order to apply the reduction,  2 files are needed:
    ogstm_South_West_I_indexes.txt
    ogstm_South_West_J_indexes.txt
    ''', formatter_class=argparse.RawTextHelpFormatter)
 
 
    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'Path of INGV forcings directory')
    parser.add_argument(   '--outputdir','-o',
                                type = str,
                                required = True,
                                help = 'Path of OGSTM forcings directory')
    parser.add_argument(   '--ingvmask','-M',
                                type = str,
                                required = True,
                                help = 'Path of the INGV mask')
 
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')
    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = False,
                                default = '19500101',
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = False,
                                default = '21000101',
                                help = 'end date in yyyymmdd format')

 
 
    return parser.parse_args()
 
args = argument()

from reducer import mesh
import numpy as np
from commons.mask import Mask
from commons.dataextractor import DataExtractor
import forcingswriter
from commons.Timelist import TimeInterval, TimeList
from commons.utils import addsep
import os

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False



def weighted_mean(values,weights):
    return (values*weights).sum()/weights.sum()

def interp2d(M2dfine, Maskfine,MaskCoarse,I,J):
    '''
    Arguments:
    * M2dfine    * 3d numpy array
    * Maskfine   * Mask object
    * MaskCoarse * Mask object
    * I          * numpy array 1d of integers
    * J          * idem
    '''
    step = I[3]-I[2]
    _, jpj, jpi = MaskCoarse.shape
        
    M2dCoarse=np.ones((jpj,jpi),np.float32)*1.e+20
    for ji in range(jpi):
        for jj in range(jpj):
            if MaskCoarse.mask[0,jj,ji]:
                i=I[ji]
                j=J[jj]
                mask  = Maskfine.mask[0,j:j+step,i:i+step]
                values=         M2dfine[j:j+step,i:i+step]
                area  =   Maskfine.area[j:j+step,i:i+step]
                M2dCoarse[jj,ji] = weighted_mean(values[mask], area[mask])
    return M2dCoarse


def interp3d(M3dfine, Maskfine,MaskCoarse,I,J):
    '''
    Degradates a 3d field defined on cell center, using area weighted mean.

    Arguments:
    * M3dfine    * 3d numpy array
    * Maskfine   * Mask object
    * MaskCoarse * Mask object
    * I          * numpy array 1d of integers, coarse_lon[k] = fine_lon[I[k]]
    * J          * idem

    Returns:
    * M3dCoarse * 3d numpy array
    '''
    step = I[3]-I[2]
    jpk, jpj, jpi = MaskCoarse.shape
        
    M3dCoarse=np.zeros((jpk,jpj,jpi),np.float32)
    for jk in range(jpk):
        for ji in range(jpi):
            for jj in range(jpj):
                if MaskCoarse.mask[jk,jj,ji]:
                    i=I[ji]
                    j=J[jj]
                    mask  = Maskfine.mask[jk,j:j+step,i:i+step]
                    values=       M3dfine[jk,j:j+step,i:i+step]
                    area  = Maskfine.area[   j:j+step,i:i+step]
                    M3dCoarse[jk,jj,ji] = weighted_mean(values[mask], area[mask])
    return M3dCoarse

def interp3d_uv(Ufine,Vfine, Maskfine,MaskCoarse,I,J):
    '''
    Generates degradated u,v fields by conserving divercence.

    Arguments:
    * Ufine      * 3d numpy array
    * Vfine      * 3d numpy array
    * Maskfine   * mesh object
    * MaskCoarse * mesh object
    * I          * numpy array 1d of integers, coarse_lon[k] = fine_lon[I[k]]
    * J          * idem

    Returns:
    * Ucoarse * 3d numpy array on coarse mask
    * Vcoarse *
    '''

    Ufine[Maskfine.umask==0] = 0.0
    Vfine[Maskfine.vmask==0] = 0.0
    step = I[3]-I[2]
    jpi = MaskCoarse.jpi
    jpj = MaskCoarse.jpj
    jpk = MaskCoarse.jpk
    UCoarse=np.zeros((jpk,jpj,jpi),np.float32)
    VCoarse=np.zeros((jpk,jpj,jpi),np.float32)

    for jk in range(jpk):
        for ji in range(jpi):
            for jj in range(jpj):               
                if (MaskCoarse.umask[jk,jj,ji]==1):
                    i=I[ji]
                    j=J[jj]
                    u =         Ufine[jk,j:j+step,i+step-1]
                    e2u= Maskfine.e2u[   j:j+step,i+step-1]
                    e3u= Maskfine.e3u_0[jk,j:j+step,i+step-1]
                    flux = (u*e2u*e3u).sum()
                    norm = MaskCoarse.e2u[jj,ji]*MaskCoarse.e3u_0[jk,jj,ji]
                    UCoarse[jk,jj,ji] = flux/norm

    for jk in range(jpk):
        for ji in range(jpi):
            for jj in range(jpj):               
                if (MaskCoarse.vmask[jk,jj,ji]==1):
                    i=I[ji]
                    j=J[jj]                     
                    v =         Vfine[jk,j+step-1,i:i+step]
                    e1v= Maskfine.e1v[   j+step-1,i:i+step]
                    e3v= Maskfine.e3v_0[jk,j+step-1,i:i+step]
                    flux = (v*e1v*e3v).sum()
                    norm = MaskCoarse.e1v[jj,ji]*MaskCoarse.e3v_0[jk,jj,ji]
                    VCoarse[jk,jj,ji] = flux/norm
    return UCoarse,VCoarse
                    

def interp2d_tau(tauxFine,tauyFine, Maskfine,MaskCoarse,I,J):
    '''
    Generates degradated taux,tauy fields by conserving stress.

    Arguments:
    * tauxFine      * 2d numpy array
    * tauyFine      * idem
    * Maskfine   * mesh object
    * MaskCoarse * mesh object
    * I          * numpy array 1d of integers, coarse_lon[k] = fine_lon[I[k]]
    * J          * idem

    Returns:
    * tauxCoarse * 2d numpy array on coarse mask
    * tauxCoarse * idem
    '''

    tauxFine[Maskfine.umask[0,:,:]==0] = 0.0
    tauyFine[Maskfine.vmask[0,:,:]==0] = 0.0
    step = I[3]-I[2]
    jpi = MaskCoarse.jpi
    jpj = MaskCoarse.jpj
    tauxCoarse=np.zeros((jpj,jpi),np.float32)
    tauyCoarse=np.zeros((jpj,jpi),np.float32)

    for ji in range(jpi):
        for jj in range(jpj):               
            if (MaskCoarse.umask[0,jj,ji]==1):
                i=I[ji]
                j=J[jj]
                u =      tauxFine[j:j+step,i+step-1]
                e2u= Maskfine.e2u[j:j+step,i+step-1]
                force = (u*e2u).sum()
                tauxCoarse[jj,ji] = force/MaskCoarse.e2u[jj,ji]

    for ji in range(jpi):
        for jj in range(jpj):               
            if (MaskCoarse.vmask[0,jj,ji]==1):
                i=I[ji]
                j=J[jj]                     
                v =      tauyFine[j+step-1,i:i+step]
                e1v= Maskfine.e1v[j+step-1,i:i+step]
                force = (v*e1v).sum()
                tauyCoarse[jj,ji] = force/MaskCoarse.e1v[jj,ji]
    return tauxCoarse,tauyCoarse


INPUTDIR = addsep(args.inputdir)
OUTDIR   = addsep(args.outputdir)

Mfine  =mesh(filename=args.ingvmask,isFreeSurface=True)
Mcoarse=mesh(filename=args.maskfile,isFreeSurface=True)
MaskINGV   = Mask(args.ingvmask)
MaskCoarse = Mask(args.maskfile)



I = np.loadtxt('South_West_I_indexes.txt',np.int32)
J = np.loadtxt('South_West_J_indexes.txt',np.int32)


TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR, "*U.nc", prefix="", dateformat="%Y%m%d")



for filename in TL.filelist[rank::nranks]:
    basename = os.path.basename(filename).rsplit("U.nc")[0]

    infileU = INPUTDIR + basename + "U.nc"
    infileV = INPUTDIR + basename + "V.nc"
    infileW = INPUTDIR + basename + "W.nc"
    infileT = INPUTDIR + basename + "T.nc"
    # same names
    outfileU = OUTDIR + basename + "U.nc"
    outfileV = OUTDIR + basename + "V.nc"
    outfileW = OUTDIR + basename + "W.nc"
    outfileT = OUTDIR + basename + "T.nc"
 
    print outfileU

    U    = DataExtractor(MaskINGV,infileU,'vozocrtx').values
    V    = DataExtractor(MaskINGV,infileV,'vomecrty').values
    TAUX = DataExtractor(MaskINGV,infileU,'sozotaux', dimvar=2).values
    TAUY = DataExtractor(MaskINGV,infileV,'sometauy', dimvar=2).values
    u,v = interp3d_uv(U, V, Mfine, Mcoarse, I, J)
    taux, tauy = interp2d_tau(TAUX, TAUY,  Mfine, Mcoarse, I, J)
    forcingswriter.writefileU(outfileU, MaskCoarse, u, taux)
    forcingswriter.writefileV(outfileV, MaskCoarse, v, tauy)



    T = DataExtractor(MaskINGV,infileT,'votemper').values
    S = DataExtractor(MaskINGV,infileT,'vosaline').values
    E = DataExtractor(MaskINGV,infileT,'sossheig', dimvar=2).values
    SW= DataExtractor(MaskINGV,infileT,'soshfldo', dimvar=2).values

    t = interp3d(T, MaskINGV, MaskCoarse, I, J)
    s = interp3d(S, MaskINGV, MaskCoarse, I, J)
    e = interp2d(E, MaskINGV, MaskCoarse, I, J)
    sw= interp2d(SW,MaskINGV, MaskCoarse, I, J)
    forcingswriter.writefileT(outfileT, MaskCoarse, t, s, e,sw)


    #W = DataExtractor(MaskINGV,infileW,'vovecrtz').values
    K = DataExtractor(MaskINGV,infileW,'votkeavt').values
    #w = interp3d(W, MaskINGV, MaskCoarse, I, J)
    k = interp3d(K, MaskINGV, MaskCoarse, I, J)

    forcingswriter.writefileW(outfileW, MaskCoarse, None, k)


