# ReadMe for msu_tbs_simple

## Overview
This package provides a simple calculate to convert climate or NWP gridded output to MSU equivalent radiances.  Since the units for the radiance output are degrees Kelvin, these radiance are also refered as "brightness temperatures".
## Input variables needed
The method used here doesn't depend on water in the atmosphere (cloud, vapor, precipitation), so the variables needed are limited.  The input is a python dictionary containing numpy arrays.  The dictionary keys are listed below.
* levels (hPa) (level).  Ordered from low pressure to higher pressure.
* temperature (K) on fixed pressure levels (lat,lon,level)
* surface_pressure  (hPa) (lat,lon)
* skin_temperature (K) and/or 2m_temperature (lat,lon).  (skin_temperature is probably more accurate for this calculation)
* land_fraction (0.0...1.0)
* sea_ice_fraction (0.0...1.0)
## Method Used
### Goals and Caveats
The goal of Method 1 is to produce reliable brightness temperatures when detailed information about moisture in the atmospheric profile is not known.  We wish to improve on the methods based on global mean weights for fixed levels both by taking into account the atmospheric pressure at the material surface, which is a measure of the total amount of absorbing material, and by making the method applicable to temperature measurements at arbitrary levels.  

### Detailed Description

We begin by calculating temperature weighting functions for each MSU channel assuming a fixed atmospheric temperature and humidity profile using oxygen and water vapor absorptivity from Rosencranz (Rosenkranz 1993; Rosenkranz 1998).  For the results presented here, we use the 1976 U.S. Standard Atmosphere (United States Committee on Extension to the Standard Atmosphere 1976) with an assumed partial pressure of water vapor, Pv , given by
$$
P_v=P_{v0} e^{-(z/z_0)}.		
$$									(1)

Here $P_{v0}$ is the vapor pressure at the surface, z is the height above the surface, and $z_0$ is the scale height, which is set to 1500 m.  $P_{v0}$ is found from the surface temperature in the U.S. Standard Atmosphere, assuming a relative humidity of 70% at the surface.  We assume that there is no liquid or solid water present.  These assumptions are arbitrary, and could easily be changed to more accurately reflect the actual vertical distribution of water vapor in the region being sampled.  The weighting functions are computed for values of the surface pressure ranging from 1100 hPa to 500 hPa with 1hPa spacing.  Separate calculations are performed for land and ocean surfaces, and for each MSU viewing angle.  The weighting function is calculated on a 1 hPa vertical grid.  We average over 13 frequencies in the pass-band of each MSU channel to account for changes in atmospheric absorption within each band.  The weighting functions for each view angle are then combined using a weighted average corresponding to the combination of views used in the MSU/AMSU product being simulated.  In Fig. 1, we plot the resulting temperature weighting functions for the TLT lower tropospheric product for 3 representative surface pressures, 1030 hPa, 990 hPa, and 700 hPa.  Obviously, pressures as low as 700 hPa are not encountered over the ocean – we include these results to illustrate the dependence of the shape of the weighting function on the surface pressure.  The figure also includes rectangular regions that represent the weight applied to the surface temperature to account for microwave emission by the surface.  These weights are also tabulated in Table 1.  Note that the relatively small decrease in surface pressure of 30 hPa (from 1000 hPa to 970 hPa) can increase the surface weight for both the land and ocean cases by more than 10%, with a corresponding, but difficult to see, decrease in atmospheric weight.  Figure 2 shows the contribution of surface weight as a function of surface pressure for land and ocean surfaces.  Typically, the land surface contribution is nearly twice the contribution from the ocean due to the much larger emissivity of land surfaces.   
	Models and radiosondes provide temperature data at a limited number of discrete reference levels.  To convert the temperatures at the reference levels to microwave brightness brightness temperatures, we need to make an assumption about the behavior of the temperature as a function of height between levels. We assume that the temperature varies linearly with height, and that the pressure is decreasing exponentially so that within the layer between two given levels the temperature T and pressure P are given by:
$$
T=T_{lower}-γ(z-z_{lower}) 
$$
$$
P=P_{lower} e^{-(z-z_{lower})/h} 							
$$
where $T_{lower}$ and $P_{lower}$ are the values at the lower level at height zlower, γ is the lapse rate, and h is the pressure scale height.  The values for h and g are found from the level data:
$$
γ=(T_{upper}-T_{lower})/(z_{lower}-z_{upper})
$$										
and
$$
h=(z_{lower}-z_{upper})/(log(P_{upper}/P_{lower})).$$										

Solving for z as a function of P and substituting, we can express T as a function of P,
$$
T=T_{lower}  (ln(P_{upper}/P))/(ln(P_{upper}/P_{lower}))+T_{upper}  (ln(P/P_{lower}))/(ln(P_{upper}/P_{lower}))
$$

a weighted combination of $T_{lower}$ and $T_{upper}$.  To calculate the weight W(n) for the level n, we sum the contributions from each 1 hPa level between the levels n-1 and n+1.  Assuming that the pressures $P_{n-1}$,$P_n$, and $P_{n+1}$ for levels n-1, n, and n+1 are whole numbers, we find that

$$
W(n)=\sum_{i=P_{n-1}}^{P_n} w(i)\frac{ln(i/P_{n-1})}{ln(P_n/P_{n-1})}+\sum_{i=P_n}^{P_{n+1}} w(i)\frac{ln(P_{n+1}/i)}{ln(P_{n+1}/P_n)}	
$$				

where w(i) is the weight for the 1 hPa level with P = i.  For the lowest layer, between the surface and the lowest reference level above the surface, the weight (corresponding to Tlower in Eq. 4)  is added into the surface emission part of the surface weight.  Thus the final “surface weight” contains contributions for the surface emissivity, and from emission from the atmosphere near the surface. For the highest layer, between the top reference level and the top of the atmosphere, all weight is assigned to the top reference level.  In Figure 3, we plot the weights assigned to a monthly average profile corresponding to the four MSU/AMSU datasets.

To calculate the brightness temperature from a temperature profile, multiply the weight for each level by the temperature and sum:
$$
T_b = W_{surface}T_{surface}+\sum_{i=0}^{n} w(i)T(i)   + W_{space}*2.73
$$


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


