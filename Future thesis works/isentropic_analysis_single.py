# Program for implementing isentropic analysis (for single time only)

import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, time):
    dataset = xr.open_dataset(nc_file)
    
    if variable_name == 'VAR_100U' or variable_name == 'VAR_100V':
        extracted_data = dataset[variable_name].sel(time=time)
        dataset.close()
    else:
        extracted_data = dataset[variable_name].sel(time=time, level=slice(100, 1000))
        dataset.close()

    return extracted_data

# FUNCTION 2: Computing wind speed (magnitude of the wind velocity vectors), converted to km/h
def wind_speed_calc(u, v):
    return np.sqrt((u*3.6)**2 + (v*3.6)**2)

# FUNCTION 3: Converting RGB color values into values usable by Matplotlib as a colormap
def read_rgb_file(rgb_file_path):
    with open(rgb_file_path, 'r') as open_file:
        rgb_lines = open_file.readlines()

    converted_colors = []
    for rgb_line in rgb_lines:
        # Split the line into components
        rgb_components = rgb_line.split()

        # Convert the RGB values to the range 0-1 and add to the list
        r, g, b = [int(c) / 255 for c in rgb_components[1:4]]
        converted_colors.append((r, g, b))

    return converted_colors

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Preparing the colormap for surface wind speed
windspeed_rgb_file_path = 'C:/Users/Brian/Desktop/Future thesis works/windspeed_color.rgb'
windspeed_colors = read_rgb_file(windspeed_rgb_file_path)
windspeed_colormap = mcolors.ListedColormap(windspeed_colors) # Colormap for wind speed

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

# isentropic_level_list = [290.0, 296.0, 300.0, 302.0, 304.0, 306.0] * units.kelvin
# Best choice for tropics: 306.0 hPa

time = '2022-12-23T00:00:00'
level = 306.0 * units.kelvin

# Data extraction
temp_data = extract_data(temp_path, 'T', time)
uwind_data = extract_data(uwind_path, 'U', time)
vwind_data = extract_data(vwind_path, 'V', time)
geopot_data = extract_data(geopot_path, 'Z', time)
spechum_data = extract_data(spechum_path, 'Q', time)
omega_data = extract_data(omega_path, 'W', time)
div_data = extract_data(div_path, 'D', time)
sfc_uwind_data = extract_data(sfc_uwind_path, 'VAR_100U', time)
sfc_vwind_data = extract_data(sfc_vwind_path, 'VAR_100V', time)

# Deriving geopotential height values from the geopotential data
geopot_hgt_data = geopot_data / 9.80665

# Combining all variable data into a single DataArray
data = xr.Dataset({
    'temp': temp_data,
    'uwind': uwind_data,
    'vwind': vwind_data,
    'geopot_hgt': geopot_hgt_data,
    'spechum': spechum_data,
    'omega': omega_data,
    'div': div_data
    })

data = data.metpy.parse_cf().squeeze()

isentropic_data = mpcalc.isentropic_interpolation_as_dataset(
    level,
    data['temp'],
    data['uwind'],
    data['vwind'],
    data['spechum'],
    data['geopot_hgt'],
    data['omega'],
    data['div']
)

isentropic_data['relhum'] = mpcalc.relative_humidity_from_specific_humidity(
    isentropic_data['pressure'],
    isentropic_data['temperature'],
    isentropic_data['spechum']
).metpy.convert_units('percent')

isentropic_data['uwind'] = isentropic_data['uwind'].metpy.convert_units('kt')
isentropic_data['vwind'] = isentropic_data['vwind'].metpy.convert_units('kt')
isentropic_data['spechum'] = (isentropic_data['spechum']) * 1000

crs = ccrs.PlateCarree()
lat = isentropic_data['pressure'].metpy.latitude
lon = isentropic_data['pressure'].metpy.longitude

bounds = [(100, 160, 0, 40)]

fig = plt.figure(figsize=(17., 12.), dpi=250)
ax = fig.add_subplot(1, 1, 1, projection=crs)
ax.set_extent(*bounds, crs=ccrs.PlateCarree())

"""
spechum_contourf = ax.contourf(lon, lat, isentropic_data['spechum'].sel(isentropic_level=level),
                               levels=np.arange(0, 20, 0.001), cmap='terrain_r', transform=ccrs.PlateCarree())

relhum_contourf = ax.contourf(lon, lat, isentropic_data['relhum'].sel(isentropic_level=level),
                              levels=np.arange(0, 130, 0.01), cmap='gist_earth_r', transform=ccrs.PlateCarree())

omega_contourf = ax.contourf(lon, lat, isentropic_data['omega'].sel(isentropic_level=level),
                             levels=np.arange(-5, 5, 0.001), cmap='seismic_r', transform=ccrs.PlateCarree())

div_contourf = ax.contourf(lon, lat, isentropic_data['div'].sel(isentropic_level=level),
                             levels=np.arange(-0.00035, 0.00035, 0.000001), cmap='seismic_r', transform=ccrs.PlateCarree())
"""

gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 5), ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

ax.coastlines('10m', edgecolor='black', linewidth=2.5)

"""
isentropic_contour_levels = np.arange(0, 1000, 25)
isentropic_contour = ax.contour(lon, lat, isentropic_data['pressure'].sel(isentropic_level=level),
isentropic_contour_levels, colors='k', linewidths=1.5, linestyles='solid',
transform=ccrs.PlateCarree())
isentropic_contour.clabel(fontsize=15, inline=1, inline_spacing=7, fmt='%i', rightside_up=True,
                          use_clabeltext=True)

spechum_colorbar = fig.colorbar(spechum_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                                extendrect='False', ticks=np.arange(0, 20, 2))
spechum_colorbar.set_label('Specific Humidity (g/kg)', size='x-large')

relhum_colorbar = fig.colorbar(relhum_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                               extendrect='False', ticks=np.arange(0, 130, 5))
relhum_colorbar.set_label('Relative Humidity (%)', size='x-large')

omega_colorbar = fig.colorbar(omega_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                              extendrect='False', ticks=np.arange(-5, 5, 0.5))
omega_colorbar.set_label('Vertical Velocity (Pa/s)', size='x-large')

div_colorbar = fig.colorbar(div_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                            extendrect='False', ticks=np.arange(-0.00035, 0.00035, 0.00005))
div_colorbar.set_label('Divergence (1/s)', size='x-large')

# Plotting wind barbs
ax.barbs(lon.values, lat.values, isentropic_data['uwind'].sel(isentropic_level=level).values,
         isentropic_data['vwind'].sel(isentropic_level=level).values, length=6,
         regrid_shape=20, transform=ccrs.PlateCarree())
"""

# Plotting surface wind velocity
sfc_wind_speed = wind_speed_calc(sfc_uwind_data, sfc_vwind_data)
sfc_wind_speed_contourf = ax.contourf(lon, lat, sfc_wind_speed, transform=ccrs.PlateCarree(),
                                      levels=np.arange(0, 185, 1), cmap=windspeed_colormap)
streamlines = ax.streamplot(lon, lat, sfc_uwind_data, sfc_vwind_data, color='k', density=2.5, arrowsize=2, arrowstyle='->')

sfc_wind_speed_colorbar = fig.colorbar(sfc_wind_speed_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                               extendrect='False', ticks=np.arange(0, 185, 10))
sfc_wind_speed_colorbar.set_label('Wind Speed (km/h)', size='x-large')

# Setting the plot titles
"""
ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Specific Humidity (g/kg) | '
             f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')

ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Relative Humidity (%) | '
             f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')

ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Vertical Velocity (Pa/s) | '
             f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')

ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Divergence (1/s) | '
             f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')
"""

ax.set_title('Surface Wind Velocity (km/h) | '
             f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')

fig.tight_layout()

plt.show()

# End of program