import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates an integral on a specific subbasin near dardanelles
                                and on 2 specific layers (0-14) (14-22)
    - the same of aveScan -
    over a season averaged files, that saves in outputdir (just for check)
    Saves an npy file.

    E.g python subbasin_analysis -o /gpfs/work/OGS18_PRACE_P_0/OPEN_BOUNDARY/ogstm_boundary_conditions/24  -f AEG_integrals_24.npy

    ''')
    parser.add_argument(   '--outfile', '-f',
                                type = str,
                                required=True,
                                help = 'file.npy')

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required=True,
                                help = '/plx/userogs/ogstsf79/OPA/V4-dev/wrkdir/2/MODEL/AVE_FREQ_1/')
    parser.set_defaults(erase=False)
    return parser.parse_args()

args = argument()
from commons.mask import Mask
from commons.submask import SubMask
from commons.Timelist import TimeList
import numpy as np
from commons import netcdf4
from commons.dataextractor import DataExtractor
from basins.region import Polygon
from basins.basin import SimplePolygonalBasin
from commons import timerequestors
from commons.time_averagers import TimeAverager3D
from commons.layer import Layer
from layer_integral.mapbuilder import MapBuilder
from commons.utils import addsep
from analysis_config import maskfile, VARLIST, Seas_obj,  mydtype, INPUTDIR, TI

OUTDIR=addsep(args.outputdir)


p = Polygon([25.2191,25.6640,26.0211,25.9826,25.6640,25.0268,24.8620],\
            [39.9939,40.1326,39.9729,39.4404,39.2322,39.5675,40.0780])
local_sub  = SimplePolygonalBasin('dar', p,'Near Dardanelles')


LAYERLIST=[ Layer(0,14), Layer(14,50)]
nLayers = len(LAYERLIST)
nSeas = Seas_obj.numbers_season
TheMask=Mask(maskfile)
S = SubMask(local_sub,maskobject=TheMask)
mask = S.mask[0,:,:]



MODEL_AEG = np.zeros((nSeas,nLayers), mydtype)

for var in VARLIST:
    TL = TimeList.fromfilenames(TI, INPUTDIR, "ave*nc", filtervar=var)
    for iSeas in range(nSeas):
        req=timerequestors.Clim_season(iSeas,Seas_obj)
   
        ii,w = TL.select(req)
        filelist= [ TL.filelist[k] for k in ii ]
        M3d = TimeAverager3D(filelist, w, var, TheMask)
        outfile = "%s%s.%s.nc" % (OUTDIR, var, Seas_obj.SEASON_LIST_NAME[iSeas])
        netcdf4.write_3d_file(M3d, var, outfile, TheMask)
        De = DataExtractor(TheMask,rawdata=M3d)
        for ilayer, layer in enumerate(LAYERLIST): 
            integrated = MapBuilder.get_layer_average(De, layer)
            MODEL_AEG[var][iSeas,ilayer] = np.nanmean(integrated[mask])
            print "%s %s %s %f " %(Seas_obj.SEASON_LIST_NAME[iSeas], var, layer, np.nanmean(integrated[mask]) )

np.save(args.outfile, MODEL_AEG)

