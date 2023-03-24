from xml.dom import minidom
import numpy as np

xmldoc = minidom.parse("rivers.xml")

mydtype=[('id',int),('name','U100'), ('mouth', 'U100'),
         ('I',int),('J',int),('SAL', np.float32),
         ('N3n',np.float32),('N1p',np.float32),
         ('DIC',np.float32),('ALK',np.float32),
         ('POC',np.float32),('DOC',np.float32)]
RIVERS=np.zeros((52,),dtype=mydtype)

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
