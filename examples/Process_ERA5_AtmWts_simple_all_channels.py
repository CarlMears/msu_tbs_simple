import numpy as np
from pathlib import Path
from msu_tbs_simple.msu_tbs_simple import MSUTbsSimple
from era5 import read_era5_data_monthly_simple, era5_monthly_files_simple
import xarray as xr
import matplotlib.pyplot as plt

try:
    from rss_plotting.global_map import plot_global_map
    do_global_map = True
except ImportError:
    print("rss_plotting module not found. Fancy Global map plotting will be unavailable.")
    print("The RSS plotting module is available here: https://github.com/CarlMears/RSS_plotting")
    do_global_map = False

if __name__ == "__main__":
    use_skin_temperature = True # Set to False to use 2m temperature instead of skin temperature in the Tbs calculation.  This is just for testing, the skin temperature should be used for the most accurate results.
    sat = 'msu'
    channel_list = ['TLT','TMT','TTS','TLS']

    era5_read_path = Path(__file__).resolve().parent / 'example_era5_data'
    output_path = Path(__file__).resolve().parent / 'output'
    output_path.mkdir(parents=True,exist_ok=True)

    for channel in channel_list:
        msu_tbs = MSUTbsSimple(channel=channel)

        year_to_do = 2024
        month_to_do = 1
        print(f'Processing Year: {year_to_do} Month: {month_to_do} Channel: {channel}')   
        
        input_files = era5_monthly_files_simple(year_to_do, month_to_do, era5_read_path)
        model_data = read_era5_data_monthly_simple(input_files)

        tbs = msu_tbs.compute_tbs(model_data=model_data, verbose=True)

        if do_global_map:
            fig,ax = plot_global_map(np.flip(tbs,0),
                            vmin=200,
                            vmax=285,
                            cmap='plasma',
                            title=f'ERA5 {channel} Tbs, {year_to_do}-{month_to_do:02d}',
                            plt_colorbar=True)
        else:
            plt.figure(figsize=(12,6))
            plt.imshow(tbs, vmin=200, vmax=285, cmap='plasma')
            plt.colorbar(label='Brightness Temperature (K)')
            plt.title(f'ERA5 {channel} Tbs, {year_to_do}-{month_to_do:02d}')
            plt.xlabel('Longitude Index')
            plt.ylabel('Latitude Index')

        ds = xr.Dataset(
            {
                'tbs': (('lat', 'lon'), tbs, {'units': 'K', 'long_name': f'ERA5 {channel} Brightness Temperature'})
            },
            coords={
                'lat': model_data['lats'],
                'lon': model_data['lons']
            },
            attrs={
                'description': f'ERA5 {channel} Tbs computed using MSUTbsSimple',
                'source': 'ERA5 reanalysis data',
                'units': 'K'
            }
        )
        output_file = output_path / f'ERA5_{channel}_Tbs_{year_to_do}_{month_to_do:02d}.nc'
        ds.to_netcdf(output_file)
        print('Wrote file: ',output_file)
    plt.show()
        



