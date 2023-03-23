from xml.dom import minidom

xmldoc = minidom.parse("rivers.xml")

for rn in xmldoc.getElementsByTagName('river'):
    print(rn.getAttribute('name'))
    salinity = float(rn.getAttribute('salinity'))
    cn = rn.getElementsByTagName('concentrations')[0]
    for vn in cn.getElementsByTagName('var'):
        print(vn.getAttribute('name'), vn.getAttribute('value'))
    