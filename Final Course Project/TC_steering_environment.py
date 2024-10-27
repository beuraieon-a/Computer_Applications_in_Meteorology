# "DEEP-LAYER MEAN WIND: Environmental Steering of Tropical Cyclones"
# This is a program for estimating the environmental sheering of a tropical cyclone (TC)
# based on the deep-layer mean (DLM) wind, i.e. the mean wind velocity within a specific
# layer of the tropopause, the thickness of which increases as a TC intensifies. This is
# an attempt to replicate the CIMSS' DLM wind product, which was anchored on the studies
# of Velden (1993) and Velden & Leslie (1991).

import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
from geopy.distance import geodesic
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import time as t

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, level_bottom, level_top, time, dp):
    dataset = xr.open_dataset(nc_file)
    wind_data_list = []
    
    while level_bottom >= level_top:
        extracted_data = dataset[variable_name].sel(level=level_bottom, time=time)
        wind_data_list.append(extracted_data)
        level_bottom = level_bottom - dp
    
    wind_data = xr.concat(wind_data_list, dim='level')
    
    dataset.close()
    
    return wind_data

# FUNCTION 2: Calculating the mean wind of each sublayer in a given steering layer
# Note: Each sublayer mean wind is already multiplied by the value of dp for convenience
# in calculating the deep-layer mean (DLM) wind. The output of this function is
# the numerator in the calculation formula for DLM wind.
def sublayer_means(isobaric_wind, level_bottom, level_top, dp):
    sublayer_means_list = []
    
    while level_bottom > level_top:
        sublayer = f'{level_bottom}–{level_bottom - dp}'
        sublayer_wind = xr.concat([isobaric_wind.sel(level = level_bottom),
                                   isobaric_wind.sel(level = level_bottom - dp)],
                                  dim='level')
        
        # Multiplying each sublayer mean by dp
        sublayer_mean_wind = sublayer_wind.mean(dim='level') * dp
        
        sublayer_mean_wind = sublayer_mean_wind.expand_dims({'sublayer': [sublayer]})
        sublayer_means_list.append(sublayer_mean_wind)
        
        level_bottom = level_bottom - dp
    
    sublayer_means_ALL = xr.concat(sublayer_means_list, dim='sublayer')
    
    return sublayer_means_ALL

# FUNCTION 3: Computing wind speed (magnitude of the wind velocity vectors), converted to km/h
def steering_wind_speed_calc(u, v):
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

# FUNCTION 5: Reading color values from a non-RGB colormap expressed in the range [1,0]
def read_colormap_file(colormap_file_path):
    with open(colormap_file_path, 'r') as open_file:
        color_lines = open_file.readlines()

    colors = []
    for color_line in color_lines:
        # Split the line into components
        r, g, b = color_line.split()
        colors.append((r, g, b))
    
    colors.reverse()
    
    return colors

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

start_time = t.time() # Recording the starting time of the program run

# Opening files: Colormaps

windspeed_rgb_file_path = 'D:/CSci135 (ComSci in Meteorology)/steeringwindspeed_colors.rgb'
# This source color file is still in RGB format and needs conversion into acceptable
# values of within the range [0,1] that is usable by Matplotlib.

vorticity_color_file_path = 'D:/CSci135 (ComSci in Meteorology)/relativevorticity_colors.rgb'
# This source color file is already in a format usable by Matplotlib,
# i.e. within the range of [0,1].

# Preparing the colormaps
windspeed_colors = read_rgb_file(windspeed_rgb_file_path)
windspeed_colormap = mcolors.ListedColormap(windspeed_colors) # Colormap for wind speed
vorticity_colors = read_colormap_file(vorticity_color_file_path)
vorticity_colormap = mcolors.ListedColormap(vorticity_colors) # Colormap for relative vorticity

# Opening files: RSMC-Tokyo TC best track dataset & ERA5 wind velocity data (u-component and v-component)
tctrack_file = 'D:/CSci135 (ComSci in Meteorology)/saola_track.txt'
uwind_file = 'D:/CSci135 (ComSci in Meteorology)/0_ERA5 Reanalysis/ERA5_uwind.nc'
vwind_file = 'D:/CSci135 (ComSci in Meteorology)/0_ERA5 Reanalysis/ERA5_vwind.nc'
TC_name = 'SAOLA' # Type in uppercase
TC_max_category = 'Super Typhoon' # Based on DOST-PAGASA TC intensity scale, type this in full (not abbreviated)

# Extracting TC track data
tctrack_data = []

with open(tctrack_file, 'r') as file:
        for line in file.readlines():
            tctrack_data.append(str(line.strip()))

time_frame = []
tc_category = []
latitude = []
longitude = []
tc_central_pressure = []
max_sustained_winds_knots = []

for per_timestep_data in tctrack_data:
    tcdata = per_timestep_data.split()
    time_frame.append(tcdata[0])
    tc_category.append(tcdata[1])
    latitude.append(float(tcdata[2]))
    longitude.append(float(tcdata[3]))
    tc_central_pressure.append(int(tcdata[4]))
    max_sustained_winds_knots.append(int(tcdata[5]))

#################################################
#### COMPUTATION OF THE DEEP-LAYER MEAN WIND ####
#################################################
# To compute the DLM wind, the program first identifies the steering layer of
# the TC based on its central pressure (CP). The thickness of the steering layer
# increases as the central pressure drops during TC intensification. Based on the
# CIMSS DLM wind product, there are six (6) steering layers, each corresponding
# to a range of central pressure: 850–700 hPa (>= 1000 hPa CP), 850–500 hPa
# (990–999 hPa CP), 850–400 hPa (970–989 hPa CP), 850–300 hPa (950–969 hPa CP)
# 850–250 hPa (940–949 hPa CP) and 700–200 hPa (<940 hPa CP). The steering layer
# is then subdivided into sublayers that are dp = 50 hPa thick. The mass-weighted
# DLM wind is calculated by: getting the mean wind per sublayer, multiplying each
# of the sublayer mean winds by the sublayer thickness (dp = 50 hPa), getting
# their summation, and then dividing the summation by the overall thickness of
# the steering layer (DP).

bottom_level = None
top_level = None
steering_layer = None

dp = float(50)
# The sublayer change in pressure level (vertical pressure change) to be used
# in DLM wind calculation, since the pressure levels to be extracted from the
# ERA5 dataset are in constant increments of dp = 50 hPa.

count = 0
# Counter for the number of 6-hourly times within the study period,
# i.e. index of the "time_frame" list

steering_product = ['Wind velocity', 'Relative vorticity']
# These are the visualization products that this program will generate.

for product in steering_product:
    for time in time_frame:
        if(tc_central_pressure[count] >= 1000):
            bottom_level = float(850)
            top_level = float(700)
            DP = float(850-700)
            steering_layer = 'Steering layer: 850–700 hPa'
        elif(999 >= tc_central_pressure[count] >= 990):
            bottom_level = float(850)
            top_level = float(500)
            DP = float(850-500)
            steering_layer = 'Steering layer: 850–500 hPa'
        elif(989 >= tc_central_pressure[count] >= 970):
            bottom_level = float(850)
            top_level = float(400)
            DP = float(850-400)
            steering_layer = 'Steering layer: 850–400 hPa'
        elif(969 >= tc_central_pressure[count] >= 950):
            bottom_level = float(850)
            top_level = float(300)
            DP = float(850-300)
            steering_layer = 'Steering layer: 850–300 hPa'
        elif(949 >= tc_central_pressure[count] >= 940):
            bottom_level = float(850)
            top_level = float(250)
            DP = float(850-250)
            steering_layer = 'Steering layer: 850–250 hPa'
        else:
            bottom_level = float(700)
            top_level = float(200)
            DP = float(700-200)
            steering_layer = 'Steering layer: 700–200 hPa'
        
        # Extracting wind data for the pressure levels inclusive of the steering layer
        # NOTE: The extracted pressure levels from ERA5 dataset are in
        # constant increments of dp = 50 hPa.
        uwind_data = extract_data(uwind_file, 'U', bottom_level, top_level, time, dp)
        vwind_data = extract_data(vwind_file, 'V', bottom_level, top_level, time, dp)

        # Defining the steering layer and calculating the sublayer means:
        # Dividing the steering layer into sublayers with dp = 50 hPa thickness,
        # i.e. into pairs of pressure levels, then concatenating the wind data in
        # each of the included pairs of pressure levels, and getting the sublayer means.
        # NOTE: The extracted pressure levels from ERA5 dataset are in constant
        # increments of dp = 50 hPa. The sublayer means are immediately multiplied
        # by dp for convenience in calculating the deep-layer mean wind.
        sublayers_uwind = sublayer_means(uwind_data, bottom_level,
                                         top_level, dp)
        sublayers_vwind = sublayer_means(vwind_data, bottom_level,
                                         top_level, dp)
        
        # Calculating the u- and v-components of the deep-layer mean wind
        u_steeringwind = sublayers_uwind.sum(dim='sublayer') / DP
        v_steeringwind = sublayers_vwind.sum(dim='sublayer') / DP
        
        # Computing speed of the deep-layer mean wind (magnitude of the wind velocity vectors) and relative vorticity
        if(product == 'Wind velocity'):
            steering_wind_speed = steering_wind_speed_calc(u_steeringwind,
                                                           v_steeringwind) # In km/h
        else:
            rel_vorticity = mpcalc.vorticity(u_steeringwind*units('m/s'),
                                             v_steeringwind*units('m/s')) # In 1/s
            rel_vorticity = rel_vorticity / (10**(-5))
        
        # Initializing the map
        if(product == 'Wind velocity'):
            x, y = np.meshgrid(steering_wind_speed['longitude'],
                               steering_wind_speed['latitude'])
        else:
            x, y = np.meshgrid(rel_vorticity['longitude'],
                               rel_vorticity['latitude'])
            
        fig = plt.figure(figsize=(25., 25.), dpi=250)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([90, 180, -10, 50], ccrs.PlateCarree()) # Full version
        # ax.set_extent([100, 150, -5, 40], ccrs.PlateCarree()) # Philippine region
        
        # Add gridlines
        gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 10),
                                 ylocs=np.arange(-90, 91, 10), color='gray', linestyle='--') # Full version
        # gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 5),
                                 # ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--') # Philippine region
        gridlines.top_labels = False
        gridlines.right_labels = False
        gridlines.xlabel_style = {'size': 20}
        gridlines.ylabel_style = {'size': 20}
        
        # Adding the boundary of the Philippine Area of Responsibility (PAR)
        PAR_longitudes = [120, 135, 135, 115, 115, 120, 120]
        PAR_latitudes = [25, 25, 5, 5, 15, 21, 25]
        ax.plot(PAR_longitudes, PAR_latitudes, color='black', linewidth=3,
                transform=ccrs.PlateCarree(), linestyle='--')
        
        # Visualizing wind speed and magnitude of relative vorticity
        if(product == 'Wind velocity'):
            wind_speed_visuals = ax.contourf(x, y, steering_wind_speed,
                                             transform=ccrs.PlateCarree(),
                                             levels=np.arange(0, 200, 1),
                                             extend='max', cmap=windspeed_colormap)
        else:
            rel_vorticity_visuals = ax.contourf(x, y, rel_vorticity,
                                                transform=ccrs.PlateCarree(),
                                                levels=np.arange(-40, 40, 1),
                                                extend='both', cmap=vorticity_colormap)
        
        # Adding the coastlines
        ax.coastlines('10m', edgecolor='black', linewidth=2.5)
        
        # Adding streamlines for wind direction
        streamlines = ax.streamplot(x, y, u_steeringwind, v_steeringwind,
                                    color='k', density=2.5, arrowsize=3,
                                    arrowstyle='->')
        
        # Plotting the TC's current location and 6-hour displacement
        latitude_fordisplacement = [latitude[count], latitude[count+1]]
        longitude_fordisplacement = [longitude[count], longitude[count+1]]
        latitude_forcurrentposition = latitude[count]
        longitude_forcurrentposition = longitude[count]
        
        ax.plot(longitude_fordisplacement, latitude_fordisplacement, 'k-X',
                transform=ccrs.PlateCarree(), linewidth=5, markersize=15)
        ax.plot(longitude_forcurrentposition, latitude_forcurrentposition, 'ro',
                transform=ccrs.PlateCarree(), markersize=15)
        
        # Calculating the TC's average forward speed
        current_position = (latitude[count], longitude[count])
        position_6hrlater = (latitude[count+1], longitude[count+1])
        
        displacement_6hr = geodesic(current_position, position_6hrlater)
        forward_speed = round(displacement_6hr.kilometers / 6, 1)
        
        # Adding a colorbar
        if(product == 'Wind velocity'):
            cax = fig.add_axes([0.92, 0.24, 0.02, 0.5]) # Full version
            # cax = fig.add_axes([0.92, 0.15, 0.02, 0.68]) # Philippine region
            norm = Normalize(vmin=0, vmax=200)
            colorbar = ColorbarBase(cax, cmap=windspeed_colormap, norm=norm,
                                    extend='max', ticks=np.arange(0, 201, 10))
            colorbar.set_label('DLM wind speed (km/h)', fontsize=30)
            colorbar.ax.tick_params(labelsize=20)
        else:
            cax = fig.add_axes([0.92, 0.24, 0.02, 0.5]) # Full version
            # cax = fig.add_axes([0.92, 0.15, 0.02, 0.68]) # Philippine region
            norm = Normalize(vmin=-40, vmax=40)
            colorbar = ColorbarBase(cax, cmap=vorticity_colormap, norm=norm,
                                    extend='both', ticks=np.arange(-40, 41, 4))
            colorbar.set_label('Relative vorticity (1/s, ×10$^{-5}$)', fontsize=30)
            colorbar.ax.tick_params(labelsize=20)
        
        # Adding plot title and captions
        title = 'Tropical cyclone environmental steering (deep-layer mean wind)'
        caption1 = '[' + time[:10] + ' ' + time[11:16] + ' UTC, ' + steering_layer + ']'
        caption2 = f'{TC_max_category} {TC_name} @ CentPress: ' + str(tc_central_pressure[count]) + ' hPa, 10-min MSW: ' + str(round((max_sustained_winds_knots[count]*1.852)/5)*5) + ' km/h, Ave. forward speed: ' + str(forward_speed) + ' km/h'
        
        # Full version
        plt.suptitle(title, x=0.125, y=0.813, fontsize=40, ha='left', va='top')
        plt.figtext(0.125, 0.788, caption1, fontsize=30, ha='left', va='top')
        plt.figtext(0.125, 0.770, caption2, fontsize=24, ha='left', va='top')

        """
        # Philippine region
        plt.suptitle(title, x=0.125, y=0.903, fontsize=40, ha='left', va='top')
        plt.figtext(0.125, 0.878, caption1, fontsize=30, ha='left', va='top')
        plt.figtext(0.125, 0.860, caption2, fontsize=24, ha='left', va='top')
        """
        
        plt.show()
    
        # Printing the status of rendering the visualization
        if(product == 'Wind velocity'):
            print('\n', product,'>> Slide', count+1, 'DONE!')
        else:
            print('\n', product,'>> Slide', count+1, 'DONE!')
        
        count = count + 1
        
        # Should stop at index = len(time_frame) - 1, this is to make sure the program ends without incurring any errors
        if(count == (len(time_frame) - 1)):
            count = 0
            break

end_time = t.time() # Recording the ending time of the program run

# Calculating and printing the program runtime
time_elapsed = end_time - start_time
hours, rem = divmod(time_elapsed, 3600)
minutes, seconds = divmod(rem, 60)
seconds, milliseconds = divmod(seconds, 1)
milliseconds *= 1000

print(f'\n PROGRAM RUNTIME: {int(hours)}:{int(minutes)}:{int(seconds)}:{int(milliseconds)}')

# End of program
