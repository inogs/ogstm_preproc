import cdsapi

c = cdsapi.Client()

c.retrieve(
    'efas-forecast',
    {
        'format': 'netcdf4.zip',
        'system_version': 'operational',
        'originating_centre': 'ecmwf',
        'product_type': 'high_resolution_forecast',
        'variable': 'river_discharge_in_the_last_6_hours',
        'model_levels': 'surface_level',
        'year': '2020',
        'month': '11',
        'day': [
            '21', '22', '23',
        ],
        'time': [
            '00:00', '12:00',
        ],
        'leadtime_hour': '72',
        'area': [
            48, 6, 35,
            25,
        ],
    },
    'download.netcdf4.zip')
