import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates 4 seasonal images
    with subbasin errorbar, for two layers.
    Each errorbar has:
    - a line for percentiles 10-90
    - a box for percentiles 25-75
    - a star for percentile 50
    Example
    https://medeaf.inogs.it/internal-validation/Benchmark/votkeavt/synthesis/Percentiles.fall.png

    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'percentiles dir')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--runname', '-n',
                                type = str,
                                required = True,
                                help = 'benchmark')

    return parser.parse_args()

args = argument()


import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
from basins import V2 as OGS
from commons import season
from commons.utils import addsep
from commons.layer import Layer


color='b'
INPUTDIR=addsep(args.inputdir)
RunName = args.runname
OUTPUTDIR=addsep(args.outdir)


SUBlist = OGS.P.basin_list
nSub = len(SUBlist)
X=np.arange(nSub)
e=0.1
S = season.season()

perc=[1,5,10,25,50,75,90,95,99]
LAYERLIST = [Layer(0,150), Layer(150,500)]

#pl.close('all')

for iSeas in range(4):
    outfile = "%sPercentiles.%s.png" %(OUTPUTDIR,S.SEASON_LIST_NAME[iSeas])
    fig,axes=pl.subplots(2,1)
    fig.set_size_inches(15,10)


    for iax, layer in enumerate(LAYERLIST):
        inputfile= "%sPercentiles.%d.%s.npy" %(INPUTDIR, iSeas, layer.string())
        P=np.load(inputfile)
        ax=axes[iax]

        for isub in range(nSub):
            p1, p5,  p10, p25, p50, p75, p90, p95, p99 = np.log10(P[:,isub])
            x = isub
            if isub == 0: 
                ax.plot([x,x],[ p10,p90 ], color, label=RunName)
            else:
                ax.plot([x,x],[ p10,p90 ], color)
            ymin, ymax = p25, p75
            ax.plot([x-e, x+e, x+e, x-e, x-e], [ymin, ymin, ymax, ymax, ymin] , color)
            ax.plot(x,p50,'*', color=color)

        if iax ==0 : ax.legend()
        ax.grid(True)
        ax.set_xticks(X)
        subnames=[sub.name for sub in OGS.P]
        ax.set_xticklabels(subnames)
        ax.set_ylabel(layer.string() + ' log10(Kz) m2/s')

    fig.suptitle("Percentiles " + S.SEASON_LIST_NAME[iSeas])
    fig.savefig(outfile)
    pl.close(fig)


