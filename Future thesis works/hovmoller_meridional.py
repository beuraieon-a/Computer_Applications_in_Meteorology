# Program for generating a meridional Hovmöller diagram

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name):
    dataset = xr.open_dataset(nc_file)
    
    if variable_name == 'VAR_100U' or variable_name == 'VAR_100V':
        extracted_data = dataset[variable_name]
        dataset.close()
    else:
        extracted_data = dataset[variable_name].sel(level=1000)
        dataset.close()

    return extracted_data

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Setting input file paths
temp_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/temperature.nc'
uwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/u-wind.nc'
vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/v-wind.nc'
geopot_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/geopotential.nc'
spechum_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/specific_humidity.nc'
omega_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/omega_velocity.nc'
div_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/divergence.nc'
sfc_uwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Surface & vertical integrals/100m_u-wind.nc'
sfc_vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Surface & vertical integrals/100m_v-wind.nc'

parameter = 'pottemp'
level = 1000 * units.hPa

# Data extraction
temp_data = extract_data(temp_path, 'T')
# uwind_data = extract_data(uwind_path, 'U')
# vwind_data = extract_data(vwind_path, 'V')
# geopot_data = extract_data(geopot_path, 'Z')
# spechum_data = extract_data(spechum_path, 'Q')
omega_data = extract_data(omega_path, 'W')
# div_data = extract_data(div_path, 'D')
# sfc_uwind_data = extract_data(sfc_uwind_path, 'VAR_100U')
# sfc_vwind_data = extract_data(sfc_vwind_path, 'VAR_100V')

"""
pottemp_data = mpcalc.potential_temperature(
    level,
    temp_data * units.kelvin
    )

relhum_data = mpcalc.relative_humidity_from_specific_humidity(
    level,
    temp_data * units.kelvin,
    spechum_data * units('kg/kg')
    ).metpy.convert_units('percent')
"""

# Data slice
lon_slice = slice(127, 128)
lat_slice = slice(30, 0)

# Get data, selecting level and lat/lon slice
# data_for_hovmoller = pottemp_data.sel(latitude=lat_slice,longitude=lon_slice)
data_for_hovmoller = omega_data.sel(latitude=lat_slice,longitude=lon_slice)

print(data_for_hovmoller)
mini = data_for_hovmoller.min().values
maxi = data_for_hovmoller.max().values
print(f'\n\n min={mini}, max={maxi}')

# Compute weights and take weighted average over longitude dimension
weights = np.cos(np.deg2rad(data_for_hovmoller.longitude.values))
average_data = (data_for_hovmoller * weights[None, None, :]).sum(dim='longitude') / np.sum(weights)

# Get times and make array of datetime objects
vtimes = data_for_hovmoller.time.values.astype('datetime64[ms]').astype('O')

# Specify latitude values for chosen domain
lats = data_for_hovmoller.latitude.values

# Start figure
fig = plt.figure(figsize=(10, 13), dpi=250)

# Plot for Hovmoller diagram
ax = plt.axes()
ax.invert_yaxis()  # Reverse the time order to do oldest first

# Plot of chosen variable averaged over longitude and slightly smoothed
# clevs = np.arange(280, 304, 1) # For pottemp
# clevs = np.arange(-0.00035, 0.00035, 0.000001) # For divergence
clevs = np.arange(-0.8, 0.8, 0.0001)
cf = ax.contourf(lats, vtimes, mpcalc.smooth_n_point(
    average_data, 9, 2), clevs, cmap='seismic_r', extend='both')
# cs = ax.contour(lats, vtimes, mpcalc.smooth_n_point(
#     average_data, 9, 2), clevs, colors='k', linewidths=1)
cbar = plt.colorbar(cf, orientation='vertical', pad=0.04, aspect=50, extendrect=True,
                    ticks=np.arange(-0.8, 0.8, 0.1))
cbar.set_label('Vertical velocity (Pa/s)')

# Make some ticks and tick labels
ax.set_xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
ax.set_xlabel('Latitude (°N)')
ax.set_yticks(vtimes[4::8])
ax.set_yticklabels(vtimes[4::8])

# Set some titles
plt.title('Hovmöller Diagram | 1000 hPa vertical velocity | {0:%Y%m%d %HZ} - {1:%Y%m%d %HZ}'.format(vtimes[0], vtimes[-1]), loc='left', fontsize=10)

plt.show()

# End of program