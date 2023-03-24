from xml.dom import minidom
import numpy as np

xmldoc = minidom.parse("rivers.xml")

mydtype=[('id',int),('name','U100'), ('mouth', 'U100'),
         ('I',int),('J',int),('SAL', np.float32),
         ('N3n',np.float32),('N1p',np.float32),
         ('DIC',np.float32),('ALK',np.float32),
         ('POC',np.float32),('DOC',np.float32)]

nPoints=52
RIVERS=np.zeros((nPoints,),dtype=mydtype)

for RiverNode in xmldoc.getElementsByTagName('river'):
    name = RiverNode.getAttribute('name')
    salinity = float(RiverNode.getAttribute('salinity'))
    MouthNode=RiverNode.getElementsByTagName('Mouths')[0]
    ConcentrationsNode = RiverNode.getElementsByTagName('concentrations')[0]
    for mn in MouthNode.getElementsByTagName('mouth'):
        ID = int(mn.getAttribute('id'))
        ind= ID -1
        RIVERS[ind]['id'] = ID
        RIVERS[ind]['I']  = int(mn.getAttribute('I'))
        RIVERS[ind]['J']  = int(mn.getAttribute('J'))
        RIVERS[ind]['SAL'] = salinity
        RIVERS[ind]['name'] = name
        RIVERS[ind]['mouth'] = mn.getAttribute('name')
        for vn in ConcentrationsNode.getElementsByTagName('var'):
            varname = vn.getAttribute('name')
            RIVERS[ind][varname] = vn.getAttribute('value')

def get_indexes_by_river(rivername):
    '''
    Searches all river mouths corresponding to a specific river name
    Arguments:
    * rivername * string, es 'Po'
    Returns:
    * good * logical array, True for points corresponding to rivername
    '''
    good=RIVERS['name']==rivername
    return good

def get_indexes_by_region(region,maskobj):
    '''
    Searches all river mouths corresponding to a specific subbasin
    Arguments:
    * region * Region object, as OGS.adr1
    Returns:
    * good * logical array, True for points corresponding to region
    '''
    good = np.zeros((nPoints),bool)
    for ir in range(nPoints):
        i = RIVERS[ir]['I'] -1
        j = RIVERS[ir]['J'] -1
        lon,lat = maskobj.convert_i_j_to_lon_lat(i, j)
        good[ir] = region.is_inside(lon,lat)
    return good

