# paths
dir_out        = "/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/BC/out/"
file_mask      = "/g100_work/OGS_test2528/V13C/QUID/SETUP/PREPROC/MASK/OGS/meshmask.nc"
file_co2       = "co2_monthly_MED_surfatm.nc"
file_nutrients = "VPnutrients_CO2.nc"
file_bmask     = dir_out + "bounmask.nc"

Perseus_river = False
if Perseus_river:
    file_river     = "Perseus-4.6_38rivers_mesh24.xlsx"



simulation_start_time = 2022
simulation_end_time   = 2024

atm_from_file = True
if atm_from_file:
    file_atm       = "atmDep_16sub_openAndCoast_winAndSum.txt"
else:
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

RST_FILES ="/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/IC/RST_2018/RST*nc"


river_data_sheet = ['KM3perYR_NOBLS',
 'DIP_KTperYR_NOBLS',
 'DIN_KTperYR_NOBLS',
 'DIS_KTperYR_NOBLS',
 'DIC_KTperYR_NOBLS',
 'ALK_GmolperYR_NOBLS',
 'O2o_GmolperYR_NOBLS',
 'DOC_KTperYR_NOBLS',
 'CDOM_KTperYR_NOBLS']


