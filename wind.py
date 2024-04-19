# Program for basic visualization of surface wind velocity using gridded meteorological data

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, time):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(time=time)
    dataset.close()

    return extracted_data

# FUNCTION 2: Data extraction for time data
def extract_time(nc_file, variable_name):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].values
    extracted_data_list = extracted_data.tolist()
    
    extracted_data_list = [np.datetime64(int(dt), 'ns') for dt in extracted_data_list]
    datetime_str_list = [str(dt) for dt in extracted_data_list]
    
    dataset.close()

    return datetime_str_list

# FUNCTION 3: Computing wind speed (magnitude of the wind velocity vectors), converted to km/h
def wind_speed_calc(u, v):
    return np.sqrt((u*3.6)**2 + (v*3.6)**2)

# FUNCTION 4: Converting RGB color values into values usable by Matplotlib as a colormap
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

# Opening files: ERA5 wind velocity data (u-component and v-component)
uwind_file = 'C:/Users/Brian/Downloads/VAR_100U.nc'
vwind_file = 'C:/Users/Brian/Downloads/VAR_100V.nc'

# Preparing wind variable names as inputs to FUNCTION 1 (extract_data)
uwind_variable_name = 'VAR_100U'
vwind_variable_name = 'VAR_100V'
time_variable_name = 'time'
# time = '2022-12-24T10:00:00'

# Data extraction
time_data = extract_time(uwind_file, time_variable_name)

for time in time_data:
    uwind = extract_data(uwind_file, uwind_variable_name, time)
    vwind = extract_data(vwind_file, vwind_variable_name, time)
    
    # Wind speed calculation
    wind_speed = wind_speed_calc(uwind, vwind) # In km/h
    
    # Initializing the map
    x, y = np.meshgrid(wind_speed['longitude'], wind_speed['latitude'])
    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([105, 160, 0, 40], ccrs.PlateCarree())
    
    # Add gridlines
    gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 5), ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
    gridlines.top_labels = False
    gridlines.right_labels = False
    
    # Visualizing wind speed
    wind_speed_visuals = ax.contourf(x, y, wind_speed, transform=ccrs.PlateCarree(), levels=np.arange(0, 185, 1), cmap=windspeed_colormap)
    
    # Adding the coastlines
    ax.coastlines('10m', edgecolor='black', linewidth=2.5)
    
    # Adding streamlines for wind direction
    streamlines = ax.streamplot(x, y, uwind, vwind, color='k', density=2.5, arrowsize=3, arrowstyle='->')
    
    # Adding a colorbar
    cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    norm = Normalize(vmin=0, vmax=185)
    colorbar = ColorbarBase(cax, cmap=windspeed_colormap, norm=norm, extend='max', ticks=np.arange(0, 186, 5))
    colorbar.set_label('Wind speed (km/h)')
    
    # Adding plot title and captions
    title = '100-m wind speed and direction | ERA5 Reanalysis'
    caption = time[:10] + ' ' + time[11:16] + ' UTC'
    plt.suptitle(title, x=0.125, y=0.830, fontsize=40, ha='left', va='top')
    plt.figtext(0.125, 0.805, caption, fontsize=30, ha='left', va='top')
    
    plt.show()
    
    print('\nWind velocity >>', time[:10], time[11:16], 'DONE!')

# End of program