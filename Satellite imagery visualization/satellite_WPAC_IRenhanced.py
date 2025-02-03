# Program for visualizing infrared (IR) satellite imagery from Himawari-8/9
# using the clean IR longwave channel (Band 13, 10.4 µm), applying a
# colormap that emphasizes extremely cold cloud-top temperatures

import os
import datetime
import numpy as np
import xarray as xr
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Extracting RGB values to convert into Matplotlib-readable format
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

# FUNCTION 2: Satellite data extraction for a specific area
def extract_satdata(file_path, band, lat_bounds, lon_bounds):
    extracted_satdata = xr.open_dataset(file_path)[band].sel(latitude=lat_bounds, longitude=lon_bounds)
    return extracted_satdata

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Preparing the colormap
EIRcolors = read_rgb_file('C:/Users/charl/Desktop/cmap/enhanced_ir.rgb')
EIRcolormap = mcolors.ListedColormap(EIRcolors)

# Setting the data extraction parameters
lat_south = -10
lat_north = 50
lon_west = 95
lon_east = 180

lat_boundary = slice(lat_north, lat_south)
lon_boundary = slice(lon_west, lon_east)

folder_path = 'C:/Users/charl/Desktop/sat_latest/'
filenames = sorted([f for f in os.listdir(folder_path) if f.endswith('.nc')])

# Setting the latitude and longitude values of the TCID region
TCID_longitudes = [110, 110, 155, 155, 110]
TCID_latitudes = [0, 27, 27, 0, 0]

# Setting the latitude and longitude values of the TCAD region
TCAD_longitudes = [114, 114, 145, 145, 114]
TCAD_latitudes = [4, 27, 27, 4, 4]

# Setting the latitude and longitude values of the PAR region
PAR_longitudes = [120, 135, 135, 115, 115, 120, 120]
PAR_latitudes = [25, 25, 5, 5, 15, 21, 25]

for filename in filenames:
    # Preparing the path to the input satellite imagery file
    sat_path = os.path.join(folder_path, filename)
    
    # Extracting satellite imagery data
    IR_cloudtop = extract_satdata(sat_path, 'tbb_13', lat_boundary, lon_boundary)  

    # Extracting the observation time
    time = str(xr.open_dataset(sat_path)['start_time'].values[0])

    #### Preparing the infrared satellite imagery data ####

    # Normalizing the IR dataset within a specific range of brightness temperatures:
    # from 188.15 K (-85Â°C) to 313.15 K (40Â°C)
    IR = (IR_cloudtop - 188.15) / (313.15 - 188.15)

    # Normalizing the IR dataset to the range 0-1
    IR = np.clip(IR, 0, 1)

    #### Visualizing the satellite imagery ####

    # Initializing the map
    x, y = np.meshgrid(IR_cloudtop['longitude'], IR_cloudtop['latitude'])
    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], ccrs.PlateCarree())
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Plotting the satellite imagery
    satellite_image = ax.imshow(IR, origin='upper',
                                extent=[lon_west, lon_east, lat_south, lat_north],
                                transform=ccrs.PlateCarree(), cmap=EIRcolormap)

    # Adding the coastlines
    ax.coastlines('10m', edgecolor='black', linewidth=1)

    # Adding gridlines
    gridlines = ax.gridlines(draw_labels=False, xlocs=np.arange(-180, 181, 5),
                             ylocs=np.arange(-90, 91, 5), color='gray',
                             linestyle='--')
    xticks = np.arange(lon_west, lon_east + 1, 5)
    yticks = np.arange(lat_south, lat_north + 1, 5)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    # Plotting the TCID boundary
    TCID, = ax.plot(TCID_longitudes, TCID_latitudes, color='yellow', linewidth=1,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Plotting the TCAD boundary
    TCAD, = ax.plot(TCAD_longitudes, TCAD_latitudes, color='orange', linewidth=1,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Plotting the PAR boundary
    PAR, = ax.plot(PAR_longitudes, PAR_latitudes, color='red', linewidth=1,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Adding a colorbar
    cax = fig.add_axes([0.914, 0.223, 0.01, 0.544])
    norm = Normalize(vmin=-85, vmax=40)
    colorbar = ColorbarBase(cax, cmap=EIRcolormap,
                            norm=norm, extend='both',
                            ticks=np.arange(-85, 41, 5))
    colorbar.set_label('Brightness temperature (°C)', fontsize=25)
    colorbar.ax.tick_params(labelsize=20)

    # Adding plot title and captions
    time_truncated = time[:16]
    time_obj = datetime.datetime.strptime(time_truncated, '%Y-%m-%dT%H:%M')
    formatted_date = time_obj.strftime('%d %B %Y')
    formatted_time = time_obj.strftime('%H:%M')

    title = 'Satellite image: Enhanced infrared (Band 13, 10.4 µm) | Himawari-9'
    caption = f'{formatted_time} UTC, {formatted_date}'
    plt.suptitle(title, x=0.125, y=0.812, fontsize=35, ha='left', va='top')
    plt.figtext(0.125, 0.789, caption, fontsize=29, ha='left', va='top')

    # Showing the visualized satellite imagery
    plt.show()

# End of program