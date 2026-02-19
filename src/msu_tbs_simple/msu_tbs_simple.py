import importlib.resources  # need this to initialize RSS_surf_emiss

from typing import Dict
from pathlib import Path
from numpy.typing import NDArray
import numpy as np

from .AtmWts_method_1 import AtmWt


class MSUTbsSimple:

    def __init__(self, channel: str) -> None:
        allowed_channels = ['TLT','TMT','TTS','TLS']
        if channel not in allowed_channels:
            raise ValueError(f"Invalid channel: {channel}. Allowed channels are: {allowed_channels}")
        self.channel = channel
        AtmWt_dict = {}
        rtm_data_path = Path(__file__).resolve().parent.parent.parent / 'data' / 'wt_tables'
        sat = 'MSU'
        AtmWt_dict['ocean'] = AtmWt(channel = channel,surface = 'ocean',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)
        AtmWt_dict['land']  = AtmWt(channel = channel,surface = 'land',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)

        self.AtmWt_dict = AtmWt_dict

    def compute_tbs(self, 
                    model_data: Dict[str, NDArray[np.float32]],
                    verbose=True,use_skin_temperature=True) -> Dict[str, NDArray[np.float32]]:

        sea_ice_frac = model_data['sea_ice_fraction']   
        land_frac = model_data['land_fraction']

        sea_ice_frac[np.isnan(sea_ice_frac)] = 0.0
        sea_ice_frac[sea_ice_frac < -1.0] = 0.0
        not_ocean = land_frac + sea_ice_frac
        not_ocean[not_ocean > 1.0] = 1.0

        #assign AtmWt classes.
        AtmWt_MSU_Ocean = self.AtmWt_dict['ocean']
        AtmWt_MSU_Land  = self.AtmWt_dict['land']

        t = model_data['temperature']
        ps = model_data['surface_pressure']
        levels= model_data['levels']

        if use_skin_temperature:
            ts = model_data['skin_temperature']
        else:
            ts = model_data['surface_temperature']
            
        # calculate the Tbs and  weights.
        if verbose:
            print('Performing Ocean Calculation')
        tbs_ocean,level_wts_ocean,surface_wts_ocean,space_wts_ocean = \
            AtmWt_MSU_Ocean.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)
        if verbose:
            print('Performing Land Calculation')
        tbs_land,level_wts_land, surface_wts_land, space_wts_land   = \
            AtmWt_MSU_Land.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)

        # combine the land and ocean results together.
        tbs_combined = not_ocean * tbs_land + (1.0 - not_ocean)*tbs_ocean

        return tbs_combined

    
