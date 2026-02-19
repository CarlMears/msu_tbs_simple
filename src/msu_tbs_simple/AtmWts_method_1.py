import numpy as np
from netCDF4 import Dataset
from numba import jit
from pathlib import Path


@jit(nopython = True)
def AtmLevelWts_Numba(*,weighting_function,
                        surface_weight,
                        space_weight,
                        pressure,
                        surface_pressure,
                        temp_profiles,
                        ts,
                        num_lats,
                        num_lons,
                        num_levels,
                        ps,
                        levels):


    MAX_SURF_PRESSURE = 1100.0
    MIN_SURF_PRESSURE = 450.0
    level_wts = np.zeros((num_levels,num_lats,num_lons))
    surface_wts = np.zeros((num_lats,num_lons))
    space_wts = np.zeros((num_lats,num_lons))
    tbs = np.zeros((num_lats,num_lons))


    for ilat in np.arange(0, num_lats):
        for ilon in range(0, num_lons):
            ps_flt = ps[ilat, ilon]

            ps_int = int(np.round(ps_flt))
            if ps_int > MAX_SURF_PRESSURE:
                raise ValueError('Surface pressure > 1100 hPa')
            if ps_int < MIN_SURF_PRESSURE:
                raise ValueError('Surface pressure < 450 hPa')
            for ps_index in np.arange(0,600):
                if ps_int == surface_pressure[ps_index]:
                    break
            wt_func  = weighting_function[:, ps_index]
            surf_wt  = surface_weight[ps_index]
            space_wt = space_weight[ps_index]

            above_surf = np.where(levels < ps_int)
            lowest_level = above_surf[0][-1]

            p_lower = ps_int
            p_upper = levels[lowest_level]
            index_lower = np.where(p_lower >= pressure)[0][0]
            index_upper = np.where(p_upper >= pressure)[0][0]

            p_layer = p_lower - 0.5 - np.arange(0, index_upper - index_lower)
            T_lower_wt_temp = 1.0 - (np.log(p_layer / p_lower) / np.log(p_upper / p_lower))
            T_upper_wt_temp = np.log(p_layer / p_lower) / np.log(p_upper / p_lower)

            T_lower_wt = np.sum(T_lower_wt_temp * wt_func[index_lower:index_upper])
            T_upper_wt = np.sum(T_upper_wt_temp * wt_func[index_lower:index_upper])

            surf_wt = surf_wt + T_lower_wt
            wt_ref = np.zeros((num_levels))
            wt_ref[lowest_level] += T_upper_wt

            # Now step through the rest of the levels
            for i in np.arange(lowest_level, 0, -1):
                p_lower = levels[i]
                p_upper = levels[i - 1]
                index_lower = np.where(p_lower >= pressure)[0][0]
                index_upper = np.where(p_upper >= pressure)[0][0]

                if index_upper > index_lower:
                    p_layer = p_lower - 0.5 - np.arange(0, index_upper - index_lower)
                    T_lower_wt_temp = 1 - (np.log(p_layer / p_lower) / np.log(p_upper / p_lower))
                    T_upper_wt_temp = np.log(p_layer / p_lower) / np.log(p_upper / p_lower)
                    T_lower_wt = np.sum(T_lower_wt_temp * wt_func[index_lower:index_upper])
                    T_upper_wt = np.sum(T_upper_wt_temp * wt_func[index_lower:index_upper])
                else:
                    T_lower_wt = 0.0
                    T_upper_wt = 0.0
                # print(i, T_lower_wt, T_upper_wt)
                wt_ref[i - 1] = wt_ref[i - 1] + T_upper_wt
                wt_ref[i] = wt_ref[i] + T_lower_wt

            p_lower = levels[0]
            p_upper = pressure[-1]
            index_lower = np.where(p_lower >= pressure)[0][0]
            index_upper = np.where(p_upper >= pressure)[0][0]

            # now the very top levels
            if ((index_upper > index_lower) and (index_upper > 0) and (index_lower > 0)):
                p_layer = p_lower - 0.5 - np.arange(0, index_upper - index_lower)
                T_lower_wt_temp = 1 - (np.log(p_layer / p_lower) / np.log(p_upper / p_lower))
                T_upper_wt_temp = np.log(p_layer / p_lower) / np.log(p_upper / p_lower)
                t_lower_wt = np.sum(T_lower_wt_temp * wt_func[index_lower:index_upper])
                t_upper_wt = np.sum(T_upper_wt_temp * wt_func[index_lower:index_upper])
            else:
                t_lower_wt = 0.0
                t_upper_wt = 0.0

            # for this last level, we include the weight from the level all the way to
            # the highest pressure

            wt_ref[0] = wt_ref[0] + t_lower_wt
            wt_ref[0] = wt_ref[0] + t_upper_wt

            level_wts[:,ilat,ilon] = wt_ref
            surface_wts[ilat,ilon] = surf_wt
            space_wts[ilat,ilon]   = space_wt

            tbs[ilat,ilon] = ts[ilat,ilon]*surf_wt + np.sum(temp_profiles[:,ilat,ilon]*wt_ref) + space_wt*2.730

    return tbs,level_wts,surface_wts,space_wts

class AtmWt():
    '''Class for level weights for MSU/AMSU'''

    def __init__(self, channel='TLT',surface = 'ocean',sat='msu',RTM_Data_Path='',verbose=True):

        path = Path(RTM_Data_Path)
        nc_file = path / f'std_atmosphere_wt_function_{sat}_chan_{channel}_{surface}_by_surface_pressure.1100.V4.nc'
        if verbose:
            print('Reading: ' + str(nc_file))
        nc_fid = Dataset(nc_file, 'r')
        self.surface = surface
        pressure = np.array(nc_fid.variables['pressure'][:])  # extract/copy the data
        self.pressure = pressure
        surface_pressure = np.array(nc_fid.variables['surface_pressure'][:])
        surface_pressure = surface_pressure.astype(np.int32)
        self.surface_pressure = surface_pressure
        surface_weight = np.array(nc_fid.variables['surface_weight'][:])
        self.surface_weight = surface_weight

        space_weight = np.array(nc_fid.variables['space_weight'][:])
        self.space_weight = space_weight
        weighting_function = np.array(nc_fid.variables['weighting_function'][:,:])
        self.weighting_function = weighting_function

    def AtmLevelWts(self, temp_profiles, ps, ts,levels):

        '''

        :param temp_profiles: numpy array of temperature profiles to be converted to MSU equivalent  Assumed to be 3D [levels,lat,lon]
                                                                                                             Assumed to ordered low P to high P
                                                                                                             Units = K
        :param ps:            numpy array containing surface pressure to be  converted.  Assumed to be 2D [lat,lon].
                                                                                         Units = hPa
        :param ts:            numpy array containing surface temperature to be converted.  Assumed to be 2D [lat,lon].
                                                                                           Units K
        :param levels:        numpy array of level pressures for the profiles in temp proffiles.  Assumed to ordered low P to high P.
                                                                                                  Units = hPa
        :return values:
                tbs          numpy array of MSU equivalent brightness temperatures, 2D, [lat,lon]
                level_wts    numpy array of level weights, 3D, [level_index,lat,lon]
                surface_wts  numpy array of surface weights, 2D, [lat,lon]
                space_wts    numpy array of "space" weights, 2D, [lat,lon].  Multiplied by 2.73K in routine
        '''

        sz1 = temp_profiles.shape
        sz2 = ps.shape
        sz3 = levels.shape

        try:
            assert(sz1[0] == sz3[0])
            assert(sz1[1] == sz2[0])
            assert(sz1[2] == sz2[1])
        except AssertionError:
            raise ValueError('Array sizes do not match  in AtmLevelWts')

        tbs,level_wts,surface_wts,space_wts = AtmLevelWts_Numba(weighting_function = self.weighting_function,
                                          surface_weight = self.surface_weight,
                                          space_weight = self.space_weight,
                                          pressure = self.pressure,
                                          surface_pressure = self.surface_pressure,
                                          temp_profiles = temp_profiles,
                                          ts = ts,
                                          num_lats = sz1[1],
                                          num_lons = sz1[2],
                                          num_levels = sz1[0],
                                          ps = ps,
                                          levels = levels)
        return tbs,level_wts,surface_wts,space_wts




