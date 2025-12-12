# paths
dir_out        = "/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/qDEG_SETUP/BC/"
file_mask      = "/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/qDEG_SETUP/MASKS/meshmask_025_z125.nc"
file_co2       = "CMIP5_scenarios_RCP_CO2_mixing_ratio.nc"
file_nutrients = "VPnutrients_CO2.nc"
file_bmask     = dir_out + "bounmask.nc"
#file_river     = "Perseus-4.6_38rivers_mesh24.xlsx"

simulation_start_time = 1999
simulation_end_time   = 2000

# Atmosphere settings in Mmol/y
# origin of data should be added here
po4_wes = (  357. +   697)/2
po4_eas = (  379. +   957)/2
n3n_wes = (10042. + 72825)/2
n3n_eas = ( 6064. + 73621)/2


# bounmask resto settings
rdpmin = 1./24
rdpmax = 90.

variables=[[u'N1p', -8.0, -6.5],
 [u'N3n', -8.0,-6.5],
 [u'O2o', -8.0,-7.5],
 [u'N5s', -8.0,-6.5],
 [u'O3c', -8.0,-6.5],
 [u'O3h', -8.0,-6.5],
 [u'N6r', -8.0,-6.5]]

gib_season = (["0630-00:00:00"])

RST_FILES ="/leonardo_work/OGS23_PRACE_IT_0/ggalli00/OGSTM-BFM/qDEG_SETUP/RESTARTS/RST*nc"


#river_data_sheet = ['KM3perYR_NOBLS',
# 'DIP_KTperYR_NOBLS',
# 'DIN_KTperYR_NOBLS',
# 'DIS_KTperYR_NOBLS',
# 'DIC_KTperYR_NOBLS',
# 'ALK_GmolperYR_NOBLS',
# 'O2o_GmolperYR_NOBLS',
# 'DOC_KTperYR_NOBLS',
# 'CDOM_KTperYR_NOBLS']


