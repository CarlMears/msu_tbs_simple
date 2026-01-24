import numpy as np
from pathlib import Path
from AtmWts_method_1 import AtmWt
from era5_monthly import read_era5_monthly_means_3D, read_era5_monthly_means_2D

import xarray as xr

def calc_era5_amsu_tbs(*,year,
                       month,
                       channel,
                       use_t2m,
                       rtm_data_path,
                       era5_read_path,
                       resolution):
    
    if resolution == '0.5/0.625':
        shift_value = 288
    elif resolution == '2.5/2.5':
        shift_value = 72
    else:
        raise ValueError('resolution must be "0.5/0.625" or "2.5/2.5"')
    
    era5_read_path = Path(era5_read_path)
    d = read_era5_monthly_means_3D(year=year, 
                                month=month, 
                                variable='temperature',
                                era5_path=era5_read_path / '3D',
                                resolution = resolution
                                )
    
    t = (d['T'][0, :, :, :]).astype(np.float32)        #Temperature in Kelvin
    levels = (d['levels'][:]).astype(np.float32)       #Pressure Levels in hPa

    d_levels = np.diff(levels)
    if np.any(d_levels < 0.0):
        print('Warning:  levels are not in ascending order.')
        if np.all(d_levels <= 0.0):
            print('Levels are in reverse Order: Reversing levels and temperature array.')
            t = np.flip(t, axis=0)
            levels = np.flip(levels)
        else:
            raise ValueError('Levels are not monotonic.  Cannot proceed.')
        
    ps = read_era5_monthly_means_2D(year=year, 
                                    month=month,
                                    variable='surface_pressure',
                                    era5_path=era5_read_path / '2D',
                                    resolution = resolution
                                    )['PS'][0, :, :]
    if np.nanmax(ps) > 2000.0:
        print('Warning: Surface pressure values appear to be in Pa, converting to hPa.')
        ps = ps / 100.0
    ps = ps.astype(np.float32)
    
    if use_t2m:   #both types of ts are in Kelvin
        ts = read_era5_monthly_means_2D(year=year, 
                                        month=month, 
                                        variable='2m_temperature',
                                        era5_path=era5_read_path / '2D',
                                        resolution = resolution
                                        )['T2m'][0, :, :]
            #this ts is the 2m air temperature.
            #  advantages of use:  more closelt tied to observations
            #  disadvantage:  Not really what the satellite sees.  The difference between T2m and Tskin can be large
            #                 under daytime and night time clear sky conditions
    else:
        ts  = read_era5_monthly_means_2D(year=year, 
                                        month=month, 
                                        variable='skin_temperature',
                                        era5_path=era5_read_path / '2D',
                                        resolution = resolution
                                        )['TSkin'][0, :, :]
            # this ts is the skin temperature as relevant to long-wave IR
            # advantage of use:  closer to microwave skin temperature, at least for moist soils
            # disadvantage:  appears to be a sort of free parameter in the model which is adjusted to satisfy
            #                energy balance
            #                for dry soil, Longwave IR and MW penetration depth are very different.
        ts = ts.astype(np.float32)



    land_frac = read_era5_monthly_means_2D(year=year, month=month, variable = 'land_sea_mask',
                                    era5_path=era5_read_path / '2D',
                                    resolution = resolution)['land_frac'][0, :, :]

    sea_ice_frac = read_era5_monthly_means_2D(year=year, month=month, variable = 'sea_ice_cover',
                                    era5_path=era5_read_path / '2D',
                                    resolution = resolution)['SeaIce'][0, :, :]

    sea_ice_frac[np.isnan(sea_ice_frac)] = 0.0
    sea_ice_frac[sea_ice_frac < -1.0] = 0.0
    not_ocean = land_frac + sea_ice_frac
    not_ocean[not_ocean > 1.0] = 1.0

    #initialize AtmWt classes.
    AtmWt_MSU_Ocean = AtmWt_dict[channel]['ocean']
    AtmWt_MSU_Land  = AtmWt_dict[channel]['land']


    # calculate the Tbs and  weights.
    print('Performing Ocean Calculation')
    tbs_ocean,level_wts_ocean,surface_wts_ocean,space_wts_ocean = \
        AtmWt_MSU_Ocean.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)
    print('Performing Land Calculation')
    tbs_land,level_wts_land, surface_wts_land, space_wts_land   = \
        AtmWt_MSU_Land.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)

    # combine the land and ocean results together.
    tbs_combined = not_ocean * tbs_land + (1.0 - not_ocean)*tbs_ocean
    tbs_combined = np.roll(np.flipud(tbs_combined),shift=shift_value,axis=1)

    lats = d['lats'][::-1]
    lons = np.roll(d['lons'],shift=shift_value)
    lons[lons >= 180.0] = lons[lons >= 180.0] - 360.0 
    tbs_xr = xr.Dataset({'tbs':(['lat','lon'],tbs_combined)},
                        coords={'lat':lats,'lon':lons},
                        attrs={'units':'K','long_name':'Simulated MSU/AMSU Tbs',
                               'channel':channel,
                               'year':f'{year:04d}',
                               'month':f'{month:02d}',
                               'model':'ERA5'})

    return tbs_xr

if __name__ == "__main__":

    resolution = '2.5/2.5'
    use_t2m = False
    sat = 'msu'
    channel_list = ['TLT']
    # initialize the AtmWt classes for all channels and surfaces
    AtmWt_dict = {}
    rtm_data_path = Path(__file__).resolve().parent / 'data' / 'wt_tables'
    for channel in channel_list:
        AtmWt_dict[channel] = {}
        AtmWt_dict[channel]['ocean'] = AtmWt(channel = channel,surface = 'ocean',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)
        AtmWt_dict[channel]['land']  = AtmWt(channel = channel,surface = 'land',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)

    if resolution == '1.0/1.0':
        resolution_str = '360_181'
    elif resolution == '2.5/2.5':
        resolution_str = '144_73'
    else:
        raise ValueError('resolution must be "1.0/1.0" or "2.5/2.5"')

    era5_read_path = Path(__file__).resolve().parent / 'example_era5_data'
    output_path = Path(__file__).resolve().parent / 'output'
    output_path.mkdir(parents=True,exist_ok=True)

    year_to_do = 2000
    month_to_do = 1
    
    for channel in channel_list:
        print(f'Processing Year: {year_to_do} Month: {month_to_do} Channel: {channel}')   
        output_file = output_path / resolution_str / f'{channel}' / f'{year_to_do:04d}' 
        if use_t2m:
            output_file = output_file / f'era5_{channel}_{year_to_do:04d}-{month_to_do:02d}.t2m.{sat}.nc'
        else:
            output_file = output_file / f'era5_{channel}_{year_to_do:04d}-{month_to_do:02d}.skt.{sat}.nc'
        output_file.parent.mkdir(parents=True,exist_ok=True) 

        tbs = calc_era5_amsu_tbs(year=year_to_do,
                                    month=month_to_do,
                                    channel=channel,
                                    use_t2m = use_t2m,
                                    rtm_data_path=rtm_data_path,
                                    era5_read_path=era5_read_path,
                                    resolution = resolution)
            
        tbs.to_netcdf(output_file)
        print('Wrote file: ',output_file)
        



