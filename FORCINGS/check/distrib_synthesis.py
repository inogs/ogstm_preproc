from bitsea.commons.layer import Layer
from bitsea.basins import V2 as OGS
import numpy as np
from bitsea.commons.utils import writetable

LayerList=[Layer(0,50), Layer(50,100), Layer(100,150), Layer(150,200), Layer(200,300), Layer(300,400), Layer(400,500)]

seasonlist=['winter','spring','summer','fall']
nLayers=len(LayerList)

#  <  1.e-4        |     > 1.0
# win spr sum fal  | win spr sum fal
OUTDIR="/g100_work/OGS_devC/Benchmark/pub/Benchmark/votkeavt/histogram_synthesis/"
INPUTDIR="/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/FORCINGS/metrics/output/"
for isub, sub in enumerate(OGS.P):
    outfile="%s%s.txt" %(OUTDIR,sub.name)
    print(outfile)
    M=np.zeros((nLayers,8))
    for iseason, season in enumerate(seasonlist):
        for ilayer, layer in enumerate(LayerList):
            inputfile="%sHistograms.%s.%s.npy" %(INPUTDIR,season,layer.string())
            A=np.load(inputfile)
            nValues=A[:,isub].sum()
            near_background_values=A[:4,isub].sum()
            great_values = A[-1,isub]
            M[ilayer,iseason] = near_background_values/nValues
            M[ilayer,iseason+4] = great_values/nValues 


    fid = open(outfile,'w')
    fid.write("#                         <  1.e-4                 |              > 1.0\n")
    fid.write("#                win     spr     sum     fal       |    win      spr     sum     fal\n")
    for ilayer, layer in enumerate(LayerList):
        fid.write(layer.longname() + "\t")
        np.savetxt(fid,M[ilayer,:4], fmt="%5.3f\t",newline="")
        fid.write("   |   ")
        np.savetxt(fid,M[ilayer,4:], fmt="%5.3f\t",newline="")
        fid.write("\n")
    fid.close()
    
    import sys
    #sys.exit()
            
            
        