# Program for conducting cross-sectional analysis

import pandas as pd
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from metpy.interpolate import cross_section
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
        extracted_data = dataset[variable_name].sel(time=time, level=slice(500, 1000))
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

# Opening files: Colormap
windspeed_rgb_file_path = 'C:/Users/Brian/Desktop/Future thesis works/windspeed_color.rgb'
# This source color file is still in RGB format and needs conversion into acceptable
# values of within the range [0,1] that is usable by Matplotlib.

# Preparing the colormaps
windspeed_colors = read_rgb_file(windspeed_rgb_file_path)
windspeed_colormap = mcolors.ListedColormap(windspeed_colors) # Colormap for wind speed

# Setting input file paths
temp_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/temperature.nc'
spechum_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/specific_humidity.nc'
uwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/u-wind.nc'
vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/v-wind.nc'
vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/v-wind.nc'
omega_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/omega_velocity.nc'
geopot_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/geopotential.nc'
div_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/divergence.nc'
sfc_uwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Surface & vertical integrals/100m_u-wind.nc'
sfc_vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Surface & vertical integrals/100m_v-wind.nc'

time_range = pd.date_range(start='2022-12-21T00:00:00', end='2022-12-27T00:00:00', freq='h')

for time in time_range:
    time = str(time)

    # Data extraction
    temp_data = extract_data(temp_path, 'T', time)
    spechum_data = extract_data(spechum_path, 'Q', time)
    uwind_data = extract_data(uwind_path, 'U', time)
    vwind_data = extract_data(vwind_path, 'V', time)
    omega_data = extract_data(omega_path, 'W', time)
    geopot_data = extract_data(geopot_path, 'Z', time)
    div_data = extract_data(div_path, 'D', time)
    sfc_uwind_data = extract_data(sfc_uwind_path, 'VAR_100U', time)
    sfc_vwind_data = extract_data(sfc_vwind_path, 'VAR_100V', time)
    
    # Deriving geopotential height values from the geopotential data
    geopot_hgt_data = geopot_data / 9.80665
    
    # Deriving mixing ratio from specific humidity
    mixratio_data = mpcalc.mixing_ratio_from_specific_humidity(spechum_data)
    
    # Deriving the w-wind from the omega vertical velocity
    wwind_data = mpcalc.vertical_velocity(omega_data, omega_data['level'], temp_data, mixratio_data)
    
    # Calculating the meridional-vertical ageostrophic wind
    ageos_v_data, ageos_w_data = mpcalc.ageostrophic_wind(
        geopot_hgt_data * units.m,
        vwind_data,
        wwind_data
        )
    
    # Combining all variable data into a single DataArray
    data = xr.Dataset({
        'temp': temp_data,
        'spechum': spechum_data,
        'mixratio': mixratio_data,
        'uwind': uwind_data,
        'vwind': vwind_data,
        'wwind': wwind_data,
        'omega': omega_data,
        'geopot_hgt': geopot_hgt_data,
        'div': div_data,
        'ageostrophic_v': ageos_v_data,
        'ageostrophic_w': ageos_w_data
        })
    
    data = data.metpy.parse_cf().squeeze()
    
    # Defining start and end point of cross-section
    start = (0.0, 127.5)
    end = (30.0, 127.5)
    
    # Extracting the cross-section
    cross = cross_section(data, start, end).set_coords(('latitude', 'longitude'))
    
    # Calculating potential temperature
    cross['pottemp'] = mpcalc.potential_temperature(
        cross['level'],
        cross['temp']
        )
    
    # Calculating relative humidity
    cross['relhum'] = mpcalc.relative_humidity_from_specific_humidity(
        cross['level'],
        cross['temp'],
        cross['spechum']
        )
    
    # Converting the wind data from m/s to kts
    cross['vwind'] = cross['vwind'].metpy.convert_units('knots')
    cross['wwind'] = cross['wwind'].metpy.convert_units('knots')
    
    # Converting specific humidity and mixing ratio from kg/kg to g/kg
    cross['spechum'] = (cross['spechum']) * 1000
    cross['mixratio'] = (cross['mixratio']) * 1000
    
    # Defining the figure object and primary axes
    fig = plt.figure(1, figsize=(16., 9.), dpi=250)
    ax = plt.axes()
    
    """
    # Plotting relative humidity
    relhum_contour = ax.contourf(cross['latitude'], cross['level'], cross['relhum'],
                                 levels=np.arange(0, 1.05, 0.005), cmap='YlGnBu')
    relhum_colorbar = fig.colorbar(relhum_contour)
    
    # Plotting specific humidity
    spechum_contour = ax.contourf(cross['latitude'], cross['level'], cross['spechum'],
                                  levels=np.arange(0, 20, 0.001), cmap='YlGnBu')
    spechum_colorbar = fig.colorbar(spechum_contour)
    
    # Plotting mixing ratio
    mixratio_contour = ax.contourf(cross['latitude'], cross['level'], cross['mixratio'],
                                   levels=np.arange(0, 20, 0.001), cmap='YlGnBu')
    mixratio_colorbar = fig.colorbar(mixratio_contour)
    
    # Plotting divergence
    div_contour = ax.contourf(cross['latitude'], cross['level'], cross['div'],
                                   levels=np.arange(-0.00035, 0.00035, 0.000001), cmap='seismic_r')
    div_colorbar = fig.colorbar(div_contour)
    """
    
    # Plotting omega vertical velocity (i.e. vertical velocity in pressure coordinates)
    omega_contour = ax.contourf(cross['latitude'], cross['level'], cross['omega'],
                                   levels=np.arange(-5, 5, 0.001), cmap='seismic_r')
    omega_colorbar = fig.colorbar(omega_contour)
    
    # Plotting potential temperature
    pottemp_contour = ax.contour(cross['latitude'], cross['level'], cross['pottemp'],
                                 levels=np.arange(250, 450, 2), colors='dimgray', linewidths=1)
    pottemp_contour.clabel(fontsize=10, colors='dimgray', inline=1, inline_spacing=8,
                           fmt='%i', rightside_up=True, use_clabeltext=True)
    
    """
    # Plotting winds
    wind_slc_vert = list(range(0, 16, 1))
    wind_slc_horz = slice(5, 100, 5)
    ax.barbs(cross['latitude'][wind_slc_horz], cross['level'][wind_slc_vert],
             cross['vwind'][wind_slc_vert, wind_slc_horz],
             cross['wwind'][wind_slc_vert, wind_slc_horz], color='k')
    """
    
    # Plotting ageostrophic winds
    wind_slc_vert = list(range(0, 16, 1))
    wind_slc_horz = slice(5, 100, 5)
    ax.quiver(cross['latitude'][wind_slc_horz], cross['level'][wind_slc_vert],
             cross['ageostrophic_v'][wind_slc_vert, wind_slc_horz],
             cross['ageostrophic_w'][wind_slc_vert, wind_slc_horz], color='k')
    
    # Adjusting the y-axis to be logarithmic
    ax.set_yscale('symlog')
    ax.set_yticklabels(np.arange(1000, 450, -50))
    ax.set_ylim(cross['level'].max(), cross['level'].min())
    ax.set_yticks(np.arange(1000, 450, -50))
    
    ax.set_xticks(np.arange(0, 31, 1))
    
    """
    # Plotting geopotential height at 500 hPa as an inset
    data_crs = ccrs.PlateCarree()
    ax_inset = fig.add_axes([0.106, 0.630, 0.25, 0.25], projection=data_crs)
    x, y = np.meshgrid(data['longitude'], data['latitude'])
    ax_inset.contour(x, y, data['geopot_hgt'].sel(level=500),
                     levels=np.arange(5100, 6000, 60), cmap='inferno')
    """
    
    # Plotting surface wind velocity as an inset
    sfc_wind_speed = wind_speed_calc(sfc_uwind_data, sfc_vwind_data)
    data_crs = ccrs.PlateCarree()
    ax_inset = fig.add_axes([0.085, 0.630, 0.25, 0.25], projection=data_crs)
    x, y = np.meshgrid(sfc_wind_speed['longitude'], sfc_wind_speed['latitude'])
    ax_inset.set_extent([100, 160, -10, 40], ccrs.PlateCarree())
    ax_inset.contourf(x, y, sfc_wind_speed, transform=ccrs.PlateCarree(),
                      levels=np.arange(0, 185, 1), cmap=windspeed_colormap)
    ax_inset.coastlines('10m')
    ax_inset.streamplot(x, y, sfc_uwind_data, sfc_vwind_data, color='k', density=1.5, linewidth=0.3, arrowsize=1, arrowstyle='->')
    
    # Plot the path of the cross section
    ax_inset.plot([start[1], end[1]], [start[0], end[0]], 'k-', transform=ccrs.PlateCarree(), linewidth=2)
    ax_inset.plot([start[1], end[1]], [start[0], end[0]], 'ro', transform=ccrs.PlateCarree())
    
    # Set the titles and axes labels
    ax_inset.set_title('')
    ax.set_title(f'ERA5 Cross-Section \u2013 {start} to {end} \u2013 '
                 f'Valid: {cross["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n'
                 'Potential Temperature (K), Meridional-Vertical Ageostrophic Winds (vectors), Vertical Velocity '
                 '(Pa/s)\nInset: Cross-Section Path and Surface Wind Velocity')
    ax.set_ylabel('Pressure (hPa)')
    ax.set_xlabel('Latitude (degrees north)')
    # relhum_colorbar.set_label('Relative Humidity (dimensionless)')
    # spechum_colorbar.set_label('Specific Humidity (g/kg)')
    # mixratio_colorbar.set_label('Mixing Ratio (g/kg)')
    # div_colorbar.set_label('Divergence (1/s)')
    omega_colorbar.set_label('Vertical Velocity (Pa/s)')
    
    plt.show()
    
    print('\n', time,'>> DONE!\n')

# End of program