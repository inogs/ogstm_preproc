from collections import namedtuple
from pathlib import Path
from xml.dom import minidom
import numpy as np


class InvalidRiverXMLFile(Exception):
    pass


# This tuple represent a BioGeoChemical variable whose average concentration
# is described inside the river.xml file
BGCVar = namedtuple('BGCVar', ['name', 'longname', 'units'])


def read_xml_vars(xml_dataset):
    """
    Read the field "Variables" of the rivers.xml file
    """
    xml_variables = xml_dataset.getElementsByTagName('Variables')

    if len(xml_variables) == 0:
        raise InvalidRiverXMLFile('No tag "Variables" found')

    if len(xml_variables) > 1:
        raise InvalidRiverXMLFile('Multiple tag "Variables" found')

    bgc_vars = []

    for node in xml_variables[0].childNodes:
        if node.nodeType == node.TEXT_NODE:
            continue
        if node.tagName != 'var':
            raise InvalidRiverXMLFile(
                'Tag with name "{}" found, when the only tag allowed is '
                '"var"'.format(node.tagName)
            )
        try:
            var_name = \
                node.getElementsByTagName("name")[0].firstChild.nodeValue
            var_longname = \
                node.getElementsByTagName("longname")[0].firstChild.nodeValue
            var_units = \
                node.getElementsByTagName("units")[0].firstChild.nodeValue
        except Exception:
            raise InvalidRiverXMLFile(
                'Error while reading the following var node:\n{}'.format(
                    node.toxml()
                )
            )

        current_var = BGCVar(
            name=var_name,
            longname=var_longname,
            units=var_units
        )

        bgc_vars.append(current_var)

    return tuple(bgc_vars)


CURRENT_DIR = Path(__file__).absolute().parent
xmldoc = minidom.parse(
    str(CURRENT_DIR / "rivers.xml")
)
BGC_VARS = read_xml_vars(xmldoc)


class Rivers:
    # The dtype for the array where we store the rivers
    DTYPE = [
        ('id', int), ('name', 'U100'), ('mouth', 'U100'), ('I', int),
        ('J', int), ('SAL', np.float32)
    ]
    # We add also 1 field for each bgc variable
    DTYPE.extend([(bgc_var.name, np.float32) for bgc_var in BGC_VARS])

    def __init__(self):
        # Count the number of mouths (they should be 52)
        rivers_section = xmldoc.getElementsByTagName('RIVERS')[0]
        n_points = len(rivers_section.getElementsByTagName('mouth'))

        # We want to be able to associate to each variable its corresponding
        # object
        self.__bgc_var_dict = {bgc_var.name: bgc_var for bgc_var in BGC_VARS}

        # We want to remember (for each variable) if we have data or not for
        # each mouth. Therefore, we associate to each variable a boolean array
        # with one entry for each mouth
        self.__var_defined = {
            bgc_var: np.zeros((n_points,), dtype=bool) for bgc_var in BGC_VARS
        }

        # We use this vector to check if we have already visited an id.
        # Now it is false everywhere
        visited_ids = np.zeros((n_points,), dtype=bool)

        # Finally, this is the array where we store the data
        rivers = np.zeros((n_points,), dtype=self.DTYPE)

        for river_node in xmldoc.getElementsByTagName('river'):
            name = river_node.getAttribute('name')
            salinity = float(river_node.getAttribute('salinity'))
            mouth_node = river_node.getElementsByTagName('Mouths')[0]
            concentrations_node = \
                river_node.getElementsByTagName('concentrations')[0]

            for mn in mouth_node.getElementsByTagName('mouth'):
                mouth_id = int(mn.getAttribute('id'))

                if mouth_id > n_points:
                    raise IOError(
                        'Mouth with id {} (for river {}) can not be saved '
                        'because the river arrays has length {}. This happens '
                        'because the ids are not contiguous (i.e., there are '
                        'some ids between 1 and {} that are not associated to '
                        'any mouth). This code can not handle this '
                        'situation.'.format(
                            mouth_id,
                            name,
                            n_points,
                            mouth_id
                        )
                    )

                ind = mouth_id - 1
                if visited_ids[ind]:
                    raise IOError(
                        'Mouth with id {} (river {}) has been defined multiple'
                        'times'.format(mouth_id, name)
                    )

                rivers[ind]['id'] = mouth_id
                # Minus 1 to be coherent with the fortran indices
                rivers[ind]['I'] = int(mn.getAttribute('I')) - 1
                rivers[ind]['J'] = int(mn.getAttribute('J')) - 1
                rivers[ind]['SAL'] = salinity
                rivers[ind]['name'] = name
                rivers[ind]['mouth'] = mn.getAttribute('name')
                for vn in concentrations_node.getElementsByTagName('var'):
                    varname = vn.getAttribute('name')
                    if varname not in self.__bgc_var_dict:
                        raise ValueError('Unknown variable: {}'.format(varname))
                    bgc_var = self.__bgc_var_dict[varname]

                    rivers[ind][varname] = vn.getAttribute('value')

                    assert not self.__var_defined[bgc_var][ind]
                    self.__var_defined[bgc_var][ind] = True

                visited_ids[ind] = True

        if not np.all(visited_ids):
            raise IOError(
                'Error in the rivers.xml files: some ids have not been '
                'defined: {}'.format(np.where(np.logical_not(visited_ids))[0])
            )
        self.__data = rivers
        self.__data.setflags(write=False)
        for var_mask in self.__var_defined.values():
            var_mask.setflags(write=False)

    def get_variable_mask(self, bgc_var):
        if isinstance(bgc_var, str):
            if bgc_var not in self.__bgc_var_dict:
                raise ValueError('Unknown variable "{}"'.format(bgc_var))
            bgc_var = self.__bgc_var_dict[bgc_var]
        return self.__var_defined[bgc_var]

    def __getattr__(self, attr_name):
        return getattr(self.__data, attr_name)

    def __getitem__(self, item):
        return self.__data[item]

    def __repr__(self):
        return repr(self.__data)

    def __str__(self):
        return str(self.__data)

    def __array__(self):
        return self.__data


RIVERS = Rivers()


def get_indexes_by_river(rivername):
    """
    Searches all river mouths corresponding to a specific river name
    Arguments:
    * rivername * string, es 'Po'
    Returns:
    * good * logical array, True for points corresponding to rivername
    """
    good = RIVERS['name'] == rivername
    return good


def get_indexes_by_region(region, maskobj):
    """
    Searches all river mouths corresponding to a specific subbasin
    Arguments:
    * region * Region object, as OGS.adr1
    Returns:
    * good * logical array, True for points corresponding to region
    """
    n_points = RIVERS.shape[0]
    good = np.zeros((n_points,), bool)
    for ir in range(n_points):
        i = RIVERS[ir]['I']
        j = RIVERS[ir]['J']
        lon, lat = maskobj.convert_i_j_to_lon_lat(i, j)
        good[ir] = region.is_inside(lon, lat)
    return good

