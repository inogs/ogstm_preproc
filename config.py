# paths
dir_out        = "out/"
file_mask      = "/marconi_scratch/userexternal/gbolzon0/EAS3_06/wrkdir/MODEL/meshmask.nc"
file_co2       = "CMIP5_scenarios_RCP_CO2_mixing_ratio.nc"
file_nutrients = "VPnutrients_CO2.nc"
file_bmask     = dir_out + "bounmask.nc"
file_river     = "Perseus-4.6_39rivers_mesh24.xlsx"

simulation_start_time = 2016
simulation_end_time   = 2022

# Atmosphere settings in Mmol/y
# origin of data should be added here
po4_wes = (  357. +   697)/2
po4_eas = (  379. +   957)/2
n3n_wes = (10042. + 72825)/2
n3n_eas = ( 6064. + 73621)/2

end_nudging=-6.0

# bounmask resto settings
rdpmin = 1./24
rdpmax = 90.

variables=[[u'N1p', -6.5],
 [u'N3n', -6.5],
 [u'O2o', -7.5],
 [u'N5s', -7.5],
 [u'O3c', -6.1],
 [u'O3h', -6.1],
 [u'N6r', -6.1]]

gib_season = (["0215-12:00:00","0515-12:00:00",
               "0815-12:00:00","1115-12:00:00"])


RST_FILES ="/gpfs/work/OGS_dev_1/REA_16/1995/INIT/RST.1995*nc"


river_data_sheet = ['KM3perYR_NOBLS',
 'DIP_KTperYR_NOBLS',
 'DIN_KTperYR_NOBLS',
 'DIS_KTperYR_NOBLS',
 'DIC_KTperYR_NOBLS',
 'ALK_GmolperYR_NOBLS',
 'O2o_GmolperYR_NOBLS']
