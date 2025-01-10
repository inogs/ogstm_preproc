import numpy as np
from pathlib import Path
from bitsea.commons import netcdf4
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
from IC import RSTwriter, smoother

date="20220801-12:00:00"
TheMask = Mask.from_file("/g100_work/OGS_prodC/OPA/Interim-dev/wrkdir/MODEL/meshmask.nc")

MDS_dir=Path("/g100_scratch/userexternal/gbolzon0/RA_24/RESTARTS/MDS")

Int_dir=Path("/g100_scratch/userexternal/gbolzon0/RA_24/RESTARTS/Interim_2023")
# generated with /g100_work/OGS_prodC/OPA/Interim-dev/archive/202307/MODEL/RST.20230801.tar
OUTDIR=Path("/g100_scratch/userexternal/gbolzon0/RA_24/RESTARTS/")

DICT = {"chl":["P1l","P2l","P3l","P4l"], "phyc":["P1c","P2c","P3c","P4c"]}

for var_mds in DICT.keys():
    VARLIST=DICT[var_mds]
    Int=np.zeros(TheMask.shape,np.float64)
    for var in VARLIST:
        Int += DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
    MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--PFTC-MedBFM3-MED-b20240216_re-sv05.00.nc", var_mds)[0,:]
    ratio = MDS / Int[:,:,80:]
    
    for var in VARLIST:
        Int = DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
        M = Int.copy()
        M[:,:,80:] = Int[:,:,80:]*ratio
        RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)
    
    


var="O2o"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--BIOL-MedBFM3-MED-b20240216_re-sv05.00.nc", "o2")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

var="N3n"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "no3")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

var="N1p"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "po4")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

var="N4n"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--NUTR-MedBFM3-MED-b20240216_re-sv05.00.nc", "nh4")
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

var="O3c"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--CARB-MedBFM3-MED-b20240216_re-sv05.00.nc", "dissic")*12000
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

var="O3h"
Int =DataExtractor(TheMask,Int_dir/f"RST.20230801-00:00:00.{var}.nc", f"TRN{var}").values
MDS = netcdf4.readfile(MDS_dir/ "20220731_d-OGS--CARB-MedBFM3-MED-b20240216_re-sv05.00.nc", "talk")*1000
M=Int.copy()
M[:,:,80:]=MDS
RSTwriter(OUTDIR / f"RST.{date}.{var}.nc",var,M,TheMask)

