import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl
from basins import V2 as OGS
from commons import season
from commons.layer import Layer



INPUTDIR1="/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/Ved_percentiles/"
INPUTDIR2="/g100_scratch/userexternal/gbolzon0/V10C/FORCINGS_EAS8/metrics/6H/Ved_percentiles/"
INPUTDIR3="/g100_scratch/userexternal/gbolzon0/V10C/FORCINGS_EAS8/metrics/24H/Ved_percentiles/"
D1={'color':'b', 'INPUTDIR': INPUTDIR1, 'runname': 'benchmark', 'shift':0}
D2={'color':'r', 'INPUTDIR': INPUTDIR2, 'runname': 'eas8_6h', 'shift':0.20}
D3={'color':'g', 'INPUTDIR': INPUTDIR3, 'runname': 'eas8_24h', 'shift':0.40}



OUTPUTDIR="/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/synthesis/"



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
        for run in [D1, D2, D3]:
            INPUTDIR=run['INPUTDIR']
            color=run['color']
            RunName=run['runname']
            shift =run['shift']
            inputfile= "%sPercentiles.%d.%s.npy" %(INPUTDIR, iSeas, layer.string())
            P=np.load(inputfile)
            ax=axes[iax]
    
            for isub in range(nSub):
                p1, p5,  p10, p25, p50, p75, p90, p95, p99 = np.log10(P[:,isub])
                x = isub + shift
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


