from commons.time_interval import TimeInterval
from static.Nutrients_reader import NutrientsReader
from static.Carbon_reader import CarbonReader
from instruments.var_conversions import NUTRVARS
import numpy as np
from commons import timerequestors
from commons import season

N=NutrientsReader()
C=CarbonReader()


def DatasetInfo(modelvarname):
    '''
    1) Helps the user to navigate inside static dataset,
    because some variables is in Nutrients one, some other
    are in the Carbon one.
    2) Defines the official name of variables to be used in deliverables
       statistics.
    Argument:
    * modelvarname * string, like 'N1p'
    Returns:
    * var     * variable name to access dataset
    * dataset * NutrientsReader or Carbonreader object
    '''
    if modelvarname in ['N1p','N3n','O2o','N5s']:
        var =  NUTRVARS[modelvarname]
        dataset     = N
    if modelvarname in ['O3h', 'Ac', 'ALK'] :
        var = 'ALK'
        dataset = C
    if modelvarname in ['O3c', 'DIC'] :
        var ='DICric'
        dataset = C
    if modelvarname in ['pH', 'PH'] :
        var='PHt_{T-Press-ins}'
        dataset = C
    if modelvarname == 'pCO2' :
        var='pCO2'
        dataset = C
    if modelvarname not in ['N1p','N3n','O2o','N5s','O3h', 'Ac','ALK','O3c', 'DIC', 'pH',',PH', 'pCO2' ]:
        raise ValueError("variable not in static dataset ")
    return var, dataset


def get_yearsStats(modelvarname, subbasinlist, LayerList, TI):
    '''
    Returns 
    * MEAN * [nsub, nLayers] numpy array, a basic annual mean
    * STD  * [nsub, nLayers] numpy array, relative std values
    '''
    nSub    = len(subbasinlist)
    nLayers = len(LayerList)
    MEAN    = np.zeros((nSub, nLayers, 5), np.float32)*np.nan
    STD     = np.zeros((nSub, nLayers, 5), np.float32)*np.nan
    NVALS   = np.zeros((nSub, nLayers, 5), np.float32)*np.nan
    S=season.season()
    var, Dataset = DatasetInfo(modelvarname)
    print var
    print Dataset
    NPROF_SEAS =np.zeros((17,5),np.int32)
    for isub, sub in enumerate(subbasinlist):
		# annual mean
		Profilelist =Dataset.Selector(var, TI, sub)
		# check
		#print "number of profiles in ", sub
		#print len(Profilelist)
		NPROF_SEAS[isub,0]=len(Profilelist)
		Pres  =np.zeros((0,),np.float32)
		Values=np.zeros((0,),np.float32)
		for p in Profilelist: 
			pres, profile, _ = p.read(var)
			Pres   = np.concatenate((Pres,pres))
			Values = np.concatenate((Values,profile))
		for ilayer, layer in enumerate(LayerList):
			ii = (Pres>=layer.top) & (Pres<layer.bottom)
			if (ii.sum()>=1) :
				MEAN[isub, ilayer, 0] = Values[ii].mean()
				STD[isub, ilayer, 0] = Values[ii].std()
				NVALS[isub,ilayer, 0] = ii.sum()

		# seasonal means
		for iseas in range(0,4): 
			reqSeas=timerequestors.Clim_season(iseas,S)
			Profilelist_seas=Dataset.Selector(var, reqSeas, sub)
			Profile_list_sel=[p for p in Profilelist_seas if TI.contains(p.time)]
			NPROF_SEAS[isub,iseas+1]=len(Profile_list_sel)
			Pres  =np.zeros((0,),np.float32)
			Values=np.zeros((0,),np.float32)
			for p in Profile_list_sel: 
				pres, profile, _ = p.read(var)
				Pres   = np.concatenate((Pres,pres))
				Values = np.concatenate((Values,profile))
				# check
				#if isub==7:
						#print "sub ", sub
						#print "iseas: ", iseas+1
						#print "profile", profile
			for ilayer, layer in enumerate(LayerList):
				ii = (Pres>=layer.top) & (Pres<layer.bottom)
				if (ii.sum()>=1) :
						MEAN[isub, ilayer, iseas+1] = Values[ii].mean()
						STD[isub, ilayer, iseas+1] =  Values[ii].std()
						NVALS[isub,ilayer,iseas+1] = ii.sum()
    return MEAN, STD, NVALS, NPROF_SEAS





