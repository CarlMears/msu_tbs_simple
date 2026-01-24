def read_era5_monthly_means_3D(year=1979,
                               month=1,
                               variable = 'temperature',
                               era5_path  = 'A:/ERA5/monthly/temperature_3D/',
                               resolution = '1.0/1.0'):
    import numpy as np
    from netCDF4 import Dataset

    if resolution == '1.0/1.0':
        resolution_str = '360_181'
    elif resolution == '2.5/2.5':
        resolution_str = '144_73'
    else:
        raise ValueError('resolution must be "1.0/1.0" or "2.5/2.5"')

    short_names = {
        'temperature': 'T',
        'specific humidity': 'Q',
        'specific cloud liquid water content': 'CLD',
        'geopotential' :  'G'
    }

    ecmwf_names = {
        'temperature': 't',
        'specific humidity': 'q',
        'specific cloud liquid water content': 'q',
        'geopotential' : 'z'
    }

    try:
        short_name = short_names[variable]
    except KeyError:
        raise ValueError('variable not in short name dictionary')

    try:
        ecmwf_name = ecmwf_names[variable]
    except KeyError:
        raise ValueError('variable not in ecmwf name dictionary')

    nc_file = era5_path / f'{year:04d}' / f'era5_{short_name}_3D_{resolution_str}_{year:04d}_{month:02d}.nc'

    nc_fid = Dataset(nc_file, 'r')
    data = np.array(nc_fid.variables[ecmwf_name])
    lons = np.array(nc_fid.variables['longitude'])
    lats = np.array(nc_fid.variables['latitude'])
    try:
        levels = np.array(nc_fid.variables['pressure_level'])
    except KeyError:
        levels = np.array(nc_fid.variables['level'])


    d = {short_name : data,
         'lats' : lats,
         'lons' : lons,
         'levels': levels,
         'name' : short_name}

    return d

def read_era5_monthly_means_2D(*,
                               year,
                               month,
                               variable,
                               era5_path,
                               resolution = '1.0/1.0'):

    import numpy as np
    from netCDF4 import Dataset

    if resolution == '1.0/1.0':
        resolution_str = '360_181'
    elif resolution == '2.5/2.5':
        resolution_str = '144_73'
    else:
        raise ValueError('resolution must be "1.0/1.0" or "2.5/2.5"')

    short_names = {
        'surface_pressure':'PS',
        '2m_temperature': 'T2m',
        'skin_temperature':'TSkin',
        '10m_wind_speed': 'W10',
        'sea_ice_cover': 'SeaIce',
        'land_sea_mask': 'land_frac',
        'total_column_water_vapour':'TPW',
        'orography' :'OR'
    }

    ecmwf_names = {
        'surface_pressure': 'sp',
        '2m_temperature': 't2m',
        'skin_temperature': 'skt',
        '10m_wind_speed': 'si10',
        'sea_ice_cover': 'siconc',
        'land_sea_mask': 'lsm',
        'total_column_water_vapour': 'tcwv',
        'orography':'z'
    }


    try:
        short_name = short_names[variable]
    except KeyError:
        raise ValueError('variable not in short name dictionary')

    try:
        ecmwf_name = ecmwf_names[variable]
    except KeyError:
        raise ValueError('variable not in ecmwf name dictionary')

    nc_file = era5_path / f'{year:04d}' / f'era5_{short_name}_2D_{resolution_str}_{year:04d}_{month:02d}.nc'

    nc_fid = Dataset(nc_file, 'r')
    data = np.array(nc_fid.variables[ecmwf_name])
    lons = np.array(nc_fid.variables['longitude'])
    lats = np.array(nc_fid.variables['latitude'])


    d = {short_name : data,
         'lats' : lats,
         'lons' : lons,
         'name' : short_name}

    return d

