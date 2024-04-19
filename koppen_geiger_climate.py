# Progam for visualizing the Köppen–Geiger climate classification map

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, lat_range, lon_range):
    dataset = xr.open_dataset(nc_file)
    all_data = dataset[variable_name]
    
    lats = dataset['lat'].values
    lons = dataset['lon'].values

    # Find indices corresponding to the specified lat/lon range
    lat_indices = np.where((lats >= lat_range[0]) & (lats <= lat_range[1]))[0]
    lon_indices = np.where((lons >= lon_range[0]) & (lons <= lon_range[1]))[0]

    extracted_data = all_data[lat_indices, lon_indices]

    dataset.close()

    return extracted_data

# FUNCTION 2: Converting RGB color values into values usable by Matplotlib as a colormap
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

# Preparing the colormap
colormap_file = 'C:/Users/Brian/Desktop/KoppenGeigerClimate_ver2.rgb'
koppen_geiger_colors = read_rgb_file(colormap_file)
koppen_geiger_colormap = mcolors.ListedColormap(koppen_geiger_colors)

# Opening NetCDF file of the Köppen–Geiger climate classification map
koppen_geiger_file = 'C:/Users/Brian/Desktop/Koppen-Geiger world climate chart, 0.1 deg.nc'

# Preparing inputs to FUNCTION 1 (extract_data)
lat_range = (4, 21.5)
lon_range = (114, 127)

# Data extraction
koppen_geiger_data = extract_data(koppen_geiger_file, 'kg_class', lat_range, lon_range)

# Initializing the map
x, y = np.meshgrid(koppen_geiger_data['lon'], koppen_geiger_data['lat'])
fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([123, 124.7, 9, 11.6], ccrs.PlateCarree())

# Add gridlines
gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 0.5), ylocs=np.arange(-90, 91, 0.5), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

# Visualizing the Köppen–Geiger climate classification
koppen_geiger_visuals = ax.pcolormesh(x, y, koppen_geiger_data, transform=ccrs.PlateCarree(), cmap=koppen_geiger_colormap)

# Adding the coastlines
ax.coastlines('10m', edgecolor='black', linewidth=2.5)

# End of program