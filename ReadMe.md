# ReadMe for msu_tbs_simple

## Overview
This package provides a simple calculate to convert climate or NWP gridded output to MSU equivalent radiances.  Since the units for the radiance output are degrees Kelvin, these radiance are also refered as "brightness temperatures".
## Input variables needed
The method used here doesn't depend on water in the atmosphere (cloud, vapor, precipitation), so the variables needed are limited.  The input is a python dictionary containing numpy arrays.  The dictionary keys are listed below.
* levels (hPa) (level).  Ordered from low pressure to higher pressure.
* temperature (K) on fixed pressure levels (lat,lon,level)
* surface_pressure  (hPa) (lat,lon)
* skin_temperature (K) or 2m_temperature (lat,lon).  (skin_temperature is probably more accurate for this calculation)
* land_fraction (0.0...1.0)
* sea_ice_fraction (0.0...1.0)
## Method Used

## Installation
### Requirements
* Python >= 3.9
* numpy
* numba
* netcdf4
### From the provided wheel
The wheel (msu_tbs_simple-{version}-py3-none-any.whl) is provide in the /dist directory. 

cd to the dist directory and install with pip:

python -m pip install msu_tbs_simple-0.1.1-py3-none-any.whl

### Installing a editable version

cd to the top directory for the package

python -m pip install -e .


## Example
An example script (and input data) is provided in the /examples directory.


