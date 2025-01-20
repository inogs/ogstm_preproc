import numpy as np
from pathlib import Path
from bitsea.commons import netcdf4
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
from IC import RSTwriter, smoother

date_out="20220701-00:00:00"
date_in="20230701-00:00:00"
date_in_mds="20220701"
TheMask = Mask.from_file("/g100_work/OGS_prodC/OPA/Interim-dev/wrkdir/MODEL/meshmask.nc")

MDS_dir=Path("/g100_work/OGS_devC/RA_24/PREPROC/RESTARTS/MDS")

Int_dir=Path("/g100_work/OGS_devC/RA_24/PREPROC/RESTARTS/Interim_2023")
# generated with /g100_work/OGS_prodC/OPA/Interim-dev/archive/202306/MODEL/RST.20230701.tar
OUTDIR=Path("/g100_work/OGS_devC/RA_24/PREPROC/RESTARTS/generated")

DICT = {"chl":["P1l","P2l","P3l","P4l"], "phyc":["P1c","P2c","P3c","P4c"]}

for var_mds in DICT.keys():
    VARLIST=DICT[var_mds]
    Int=np.zeros(TheMask.shape,np.float64)
    for var in VARLIST:
        Int += DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
    MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--PFTC-MedBFM3-MED-b20240216_re-sv05.00.nc", var_mds)[0,:]
    ratio = MDS / Int[:,:,80:]
    
    for var in VARLIST:
        Int = DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
        M = Int.copy()
        M[:,:,80:] = Int[:,:,80:]*ratio
        RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)
    
    


var="O2o"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--BIOL-MedBFM3-MED-b20240216_re-sv05.00.nc", "o2")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

var="N3n"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "no3")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

var="N1p"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "po4")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

var="N4n"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "nh4")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

var="O3c"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--CARB-MedBFM3-MED-b20240216_re-sv05.00.nc", "dissic")*12000
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

var="O3h"
Int =DataExtractor(TheMask,Int_dir/f"RST.{date_in}.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ f"{date_in_mds}_d-OGS--CARB-MedBFM3-MED-b20240216_re-sv05.00.nc", "talk")*1000
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date_out}.{var}.nc",var,M,TheMask)

