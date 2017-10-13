# paths
dir_out        = "out/"
file_mask      = "meshmask_469.nc"
file_co2       = "CMIP5_scenarios_RCP_CO2_mixing_ratio.nc"
file_nutrients = "VPnutrients_CO2.nc"
file_bmask     = "out/bounmask.nc"
file_river     = "river_8_on_meshgrid_v2_4.xlsx"

simulation_start_time = 1999
simulation_end_time   = 2020

# Atmosphere settings in Mmol/m3/y
# origin of data should be added here
po4_wes = (  357. +   697)/2
po4_eas = (  379. +   957)/2
n3n_wes = (10042. + 72825)/2
n3n_eas = ( 6064. + 73621)/2

end_nudging=-5.25

# bounmask resto settings
rdpmin = 1./24
rdpmax = 90.

variables=[[u'N1p', -7.5],
 [u'N3n', -7.5],
 [u'O2o', -7.5],
 [u'N5s', -7.5],
 [u'O3c', -5.5],
 [u'O3h', -5.5],
 [u'N6r', -5.5]]


river_data_sheet = ['KM3perYR_NOBLS',
 'DIP_KTperYR_NOBLS',
 'DIN_KTperYR_NOBLS',
 'DIS_KTperYR_NOBLS',
 'DIC_KTperYR_NOBLS',
 'ALK_GmolperYR_NOBLS',
 'O2o_GmolperYR_NOBLS']
