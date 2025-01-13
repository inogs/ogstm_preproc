# paths
dir_out        = "/g100_work/OGS23_PRACE_IT/grosati/NECCTON/PREPROC_ogstm/ogstm_preproc/BC/out/"
file_mask      = "/g100_work/OGS23_PRACE_IT/grosati/NECCTON/PREPROC_ogstm/ogstm_preproc/BC/meshmask.nc"
file_co2       = "CMIP5_scenarios_RCP_CO2_mixing_ratio.nc"
file_nutrients = "VPnutrients_CO2Hg.nc"
file_bmask     = "/g100_scratch/userexternal/camadio0/Neccton_hindcast1999_2022/wrkdir/MODEL/bounmask.nc"
file_river     = "Perseus-4.6_40rivers_genericmesh_withDissHg.xlsx"

simulation_start_time = 2015
simulation_end_time   = 2020

# Atmosphere settings in Mmol/y
# origin of data should be added here
po4_wes = (  357. +   697)*4./2.
po4_eas = (  379. +   957)*4./2.
n3n_wes = (10042. + 72825)/2.
n3n_eas = ( 6064. + 73621)/2.
HgII_med = 0.19    #HgII deposition Mmol/y
MMHg_med = HgII_med*0.02 #MMHg deposition Mmol/y
Hg0atm_med =  1.6         # ng/g 64. + 73621)/2.


# bounmask resto settings
rdpmin = 1./24
rdpmax = 90.

variables=[[u'N1p', -8.0, -6.5],
 [u'N3n', -8.0,-6.5],
 [u'O2o', -8.0,-7.5],
 [u'N5s', -8.0,-6.5],
 [u'O3c', -8.0,-6.5],
 [u'O3h', -8.0,-6.5],
 [u'N6r', -8.0,-6.5],
 [u'Hg0', -8.0,-6.5],
 [u'Hg2', -8.0,-6.5],
 [u'MHg',-8.0,-6.5],
 [u'DHg',-8.0,-6.5]]

gib_season = (["0630-00:00:00"])

#RST_FILES ="../PREPROC/IC/RST_2018/RST*nc"
RST_FILES ="/g100_work/OGS23_PRACE_IT/grosati/NECCTON/PREPROC_ogstm/ogstm_preproc/IC/RST_const/RST*nc"


river_data_sheet = ['KM3perYR_NOBLS',
 'DIP_KTperYR_NOBLS',
 'DIN_KTperYR_NOBLS',
 'DIS_KTperYR_NOBLS',
 'DIC_KTperYR_NOBLS',
 'ALK_GmolperYR_NOBLS',
 'O2o_GmolperYR_NOBLS',
 'HgII_KTperYR_NOBLS',
 'MMHg_KTperYR_NOBLS',
 'DOC_KTperYR_NOBLS',
 'CDOM_KTperYR_NOBLS']


