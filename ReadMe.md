# ReadMe for msu_tbs_simple

## Overview
This package provides a simple calculate to convert climate or NWP gridded output to MSU equivalent radiances.  Since the units for the radiance output are degrees Kelvin, these radiance are also refered as "brightness temperatures".
## Input variables needed
The method used here doesn't depend on water in the atmosphere (cloud, vapor, precipitation), so the variables needed are limited.  The input is a python dictionary containing numpy arrays.  The dictionary keys are listed below.
* temperature (K) on fixed pressure levels (lat,lon,level)
* surface_pressure  (hPa) (lat,lon)
* skin_temperature (K) or 2m_temperature (lat,lon).  (Skin_temperature is probably more accurate for this calculation)
* land_fraction (0.0...1.0)
* sea_ice_fraction (0.0...1.0)
