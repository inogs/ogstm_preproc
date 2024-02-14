import argparse


def read_command_line_args():
    parser = argparse.ArgumentParser(
        description='Generates PO files reading runoff from forcings'
    )
    parser.add_argument(
        '--inputdir',
        '-i',
        type=str,
        required=True,
        help='/some/path/'
    )
    parser.add_argument(
        '--outdir',
        "-o",
        type=str,
        required=True,
        help='/some/path/'
    )
    parser.add_argument(
        '--cmccmaskfile', '-M',
        type=str,
        required=True,
        help='/some/path/outmask.nc'
    )
    parser.add_argument(
        '--maskfile', '-m',
        type=str,
        required=True,
        help='''NetCDF File name of the mask.''')

    return parser.parse_args()


ARGS = read_command_line_args()


from collections import namedtuple
import numpy as np
import netCDF4
from os import path
import seawater as sw

from commons.utils import addsep
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.dataextractor import DataExtractor

from river_reader import RIVERS, BGC_VARS


# Some variables must be computed by changing the name and/or the units of
# measurements of the original ones. In this dictionary, we save the
# conversions of each variable: we associate the original name to a tuple of
# conversions; each conversion is a pair that contains a new name and the
# coefficient we use to multiply the original variable. If a variable is not in
# this dictionary, it will be copied "as is" (i.e., its conversion will be the
# tuple that contains only the pair (original_name, 1.))
VARS_CONVERSIONS = {
    'ALK': (('O3h', '1. / rho'),),
    'DIC': (('O3c', '1. / rho'),),
    'POC': (('R6c', '1. / rho'),),
    'DOC': (('R3c', '1. / rho'),),
    'N3n': (('N5s', 28.0855 / 14.0067), ('N3n', 1.))
}


FILL_VALUE = 1e20


OutputVariable = namedtuple(
    'OutputVariable',
    ('name', 'lon_positions', 'lat_positions', 'values')
)


def main():
    input_dir = addsep(ARGS.inputdir)
    output_dir = addsep(ARGS.outdir)

    cmcc_mask = Mask(ARGS.cmccmaskfile)
    ogs_mask = Mask(ARGS.maskfile)
    _,jpj,jpi=ogs_mask.shape
    _,_,JPI =cmcc_mask.shape

    # Check the size of the region of CMCC mask on the left of the strait of
    # Gibraltar that is missing in the OGS model
    mask_cut = JPI - jpi

    lon_positions = RIVERS['I']
    lat_positions = RIVERS['J']

    surface_level_height = cmcc_mask.zlevels[0]
    pressure = np.ones_like(lon_positions, dtype=np.float32)
    pressure *= surface_level_height

    TL = TimeList.fromfilenames(None, input_dir, "T*.nc", prefix="T")

    for timestep, filename in enumerate(TL.filelist):
        runoff = DataExtractor(
            cmcc_mask,
            filename,
            'sorunoff',
            dimvar=2
        ).values[lat_positions, lon_positions + mask_cut]

        surface_potential_temperature = DataExtractor(
            cmcc_mask,
            filename,
            'votemper',
            dimvar=2
        ).values[lat_positions, lon_positions + mask_cut]

        surface_temperature = sw.temp(
            s=RIVERS['SAL'],
            pt=surface_potential_temperature,
            p=pressure
        )

        density = sw.dens(
            s=RIVERS['SAL'],
            t=surface_temperature,
            p=pressure
        )

        cell_areas = cmcc_mask.area[lat_positions, lon_positions + mask_cut]

        discharge = runoff * cell_areas / density  # m^3/s

        output_data = []
        for variable in BGC_VARS:
            # If there is not the name of the variable on the VARS_CONVERSIONS
            # dict, create a dummy conversion with the same name and coefficient
            # equal to 1.
            var_conversions = VARS_CONVERSIONS.get(
                variable.name,
                ((variable.name, 1.),)
            )

            var_mask = RIVERS.get_variable_mask(variable)
            raw_var_concentration = RIVERS[variable.name][var_mask]

            # Keep only the values related to the current variable
            var_discharge = discharge[var_mask]
            var_cell_areas = cell_areas[var_mask]

            for var_name, conversion_factor_raw in var_conversions:
                if isinstance(conversion_factor_raw, str):
                    rho=density[var_mask]
                    conversion_factor = np.asarray(eval(conversion_factor_raw))
                else:
                    conversion_factor = conversion_factor_raw

                var_concentration = raw_var_concentration * conversion_factor
                var_data = var_discharge * var_concentration / var_cell_areas

                current_var = OutputVariable(
                    name=var_name,
                    lon_positions=lon_positions[var_mask],
                    lat_positions=lat_positions[var_mask],
                    values=var_data
                )
                output_data.append(current_var)

        # Now we write the content of the output_data dictionary on a netcdf
        # file
        date_str = TL.Timelist[timestep].strftime("%Y%m%d-%H:%M:%S")
        outfile = path.join(
            output_dir,
            "TIN_{}.nc".format(date_str)
        )
        print(outfile,flush=True)
        with netCDF4.Dataset(outfile, 'w') as nc_file:
            nc_file.createDimension('lon', jpi)
            nc_file.createDimension('lat', jpj)

            for output_var in output_data:
                nc_var_data=np.ones((jpj,jpi),np.float32)*FILL_VALUE

                # Fill water points with -1 to help visualization
                nc_var_data[ogs_mask.mask_at_level(0)] = -1

                lon_pos = output_var.lon_positions
                lat_pos = output_var.lat_positions
                nc_var_data[(lat_pos, lon_pos)] = output_var.values

                nc_var = nc_file.createVariable(
                    'riv_{}'.format(output_var.name),
                    datatype='f4',
                    dimensions=('lat', 'lon'),
                    fill_value=FILL_VALUE,
                    zlib=True,
                    complevel=2
                )

                nc_var[:] = nc_var_data


if __name__ == '__main__':
    main()