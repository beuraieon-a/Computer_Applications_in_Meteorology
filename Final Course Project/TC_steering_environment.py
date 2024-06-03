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
def extract_data(nc_file, variable_name, pressure_level, time):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(level=pressure_level, time=time)
    dataset.close()

    return extracted_data

# FUNCTION 2: Computing wind speed (magnitude of the wind velocity vectors), converted to km/h
def steering_wind_speed_calc(u, v):
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

# FUNCTION 4: Reading color values from a non-RGB colormap expressed in the range [1,0]
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

windspeed_rgb_file_path = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/steeringwindspeed_colors.rgb'
# This source color file is still in RGB format and needs conversion into acceptable
# values of within the range [0,1] that is usable by Matplotlib.

vorticity_color_file_path = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/relativevorticity_colors.rgb'
# This source color file is already in a format usable by Matplotlib,
# i.e. within the range of [0,1].

# Preparing the colormaps
windspeed_colors = read_rgb_file(windspeed_rgb_file_path)
windspeed_colormap = mcolors.ListedColormap(windspeed_colors) # Colormap for wind speed
vorticity_colors = read_colormap_file(vorticity_color_file_path)
vorticity_colormap = mcolors.ListedColormap(vorticity_colors) # Colormap for relative vorticity

# Opening files: RSMC-Tokyo TC best track dataset & ERA5 wind velocity data (u-component and v-component)
tctrack_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/saola_track.txt'
uwind_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/0_ERA5 Reanalysis/ERA5_uwind.nc'
vwind_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/0_ERA5 Reanalysis/ERA5_vwind.nc'

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

# Preparing wind variable names as inputs to FUNCTION 1 (extract_data)
uwind_variable_name = 'U'
vwind_variable_name = 'V'

count = 0
# Counter for the number of 6-hourly times within the study period,
# i.e. index of the "time_frame" list

steering_product = ['Wind velocity', 'Relative vorticity']
# These are the visualization products that this program will generate.

dp = float(50)
# The sublayer change in pressure level (vertical pressure change) to be used
# in DLM wind calculation, since the pressure levels to be extracted from the
# ERA5 dataset are in constant increments of dp = 50 hPa.

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

for product in steering_product:
    for time in time_frame:
        # Extracting wind data for the pressure levels inclusive of the steering layer
        # NOTE: The extracted pressure levels from ERA5 dataset are in
        # constant increments of dp = 50 hPa.
        if(tc_central_pressure[count] >= 1000):
            uwind850_data = extract_data(uwind_file, uwind_variable_name, 850, time)
            vwind850_data = extract_data(vwind_file, vwind_variable_name, 850, time)
            uwind800_data = extract_data(uwind_file, uwind_variable_name, 800, time)
            vwind800_data = extract_data(vwind_file, vwind_variable_name, 800, time)
            uwind750_data = extract_data(uwind_file, uwind_variable_name, 750, time)
            vwind750_data = extract_data(vwind_file, vwind_variable_name, 750, time)
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
        elif(999 >= tc_central_pressure[count] >= 990):
            uwind850_data = extract_data(uwind_file, uwind_variable_name, 850, time)
            vwind850_data = extract_data(vwind_file, vwind_variable_name, 850, time)
            uwind800_data = extract_data(uwind_file, uwind_variable_name, 800, time)
            vwind800_data = extract_data(vwind_file, vwind_variable_name, 800, time)
            uwind750_data = extract_data(uwind_file, uwind_variable_name, 750, time)
            vwind750_data = extract_data(vwind_file, vwind_variable_name, 750, time)
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
            uwind650_data = extract_data(uwind_file, uwind_variable_name, 650, time)
            vwind650_data = extract_data(vwind_file, vwind_variable_name, 650, time)
            uwind600_data = extract_data(uwind_file, uwind_variable_name, 600, time)
            vwind600_data = extract_data(vwind_file, vwind_variable_name, 600, time)
            uwind550_data = extract_data(uwind_file, uwind_variable_name, 550, time)
            vwind550_data = extract_data(vwind_file, vwind_variable_name, 550, time)
            uwind500_data = extract_data(uwind_file, uwind_variable_name, 500, time)
            vwind500_data = extract_data(vwind_file, vwind_variable_name, 500, time)
        elif(989 >= tc_central_pressure[count] >= 970):
            uwind850_data = extract_data(uwind_file, uwind_variable_name, 850, time)
            vwind850_data = extract_data(vwind_file, vwind_variable_name, 850, time)
            uwind800_data = extract_data(uwind_file, uwind_variable_name, 800, time)
            vwind800_data = extract_data(vwind_file, vwind_variable_name, 800, time)
            uwind750_data = extract_data(uwind_file, uwind_variable_name, 750, time)
            vwind750_data = extract_data(vwind_file, vwind_variable_name, 750, time)
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
            uwind650_data = extract_data(uwind_file, uwind_variable_name, 650, time)
            vwind650_data = extract_data(vwind_file, vwind_variable_name, 650, time)
            uwind600_data = extract_data(uwind_file, uwind_variable_name, 600, time)
            vwind600_data = extract_data(vwind_file, vwind_variable_name, 600, time)
            uwind550_data = extract_data(uwind_file, uwind_variable_name, 550, time)
            vwind550_data = extract_data(vwind_file, vwind_variable_name, 550, time)
            uwind500_data = extract_data(uwind_file, uwind_variable_name, 500, time)
            vwind500_data = extract_data(vwind_file, vwind_variable_name, 500, time)
            uwind450_data = extract_data(uwind_file, uwind_variable_name, 450, time)
            vwind450_data = extract_data(vwind_file, vwind_variable_name, 450, time)
            uwind400_data = extract_data(uwind_file, uwind_variable_name, 400, time)
            vwind400_data = extract_data(vwind_file, vwind_variable_name, 400, time)
        elif(969 >= tc_central_pressure[count] >= 950):
            uwind850_data = extract_data(uwind_file, uwind_variable_name, 850, time)
            vwind850_data = extract_data(vwind_file, vwind_variable_name, 850, time)
            uwind800_data = extract_data(uwind_file, uwind_variable_name, 800, time)
            vwind800_data = extract_data(vwind_file, vwind_variable_name, 800, time)
            uwind750_data = extract_data(uwind_file, uwind_variable_name, 750, time)
            vwind750_data = extract_data(vwind_file, vwind_variable_name, 750, time)
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
            uwind650_data = extract_data(uwind_file, uwind_variable_name, 650, time)
            vwind650_data = extract_data(vwind_file, vwind_variable_name, 650, time)
            uwind600_data = extract_data(uwind_file, uwind_variable_name, 600, time)
            vwind600_data = extract_data(vwind_file, vwind_variable_name, 600, time)
            uwind550_data = extract_data(uwind_file, uwind_variable_name, 550, time)
            vwind550_data = extract_data(vwind_file, vwind_variable_name, 550, time)
            uwind500_data = extract_data(uwind_file, uwind_variable_name, 500, time)
            vwind500_data = extract_data(vwind_file, vwind_variable_name, 500, time)
            uwind450_data = extract_data(uwind_file, uwind_variable_name, 450, time)
            vwind450_data = extract_data(vwind_file, vwind_variable_name, 450, time)
            uwind400_data = extract_data(uwind_file, uwind_variable_name, 400, time)
            vwind400_data = extract_data(vwind_file, vwind_variable_name, 400, time)
            uwind350_data = extract_data(uwind_file, uwind_variable_name, 350, time)
            vwind350_data = extract_data(vwind_file, vwind_variable_name, 350, time)
            uwind300_data = extract_data(uwind_file, uwind_variable_name, 300, time)
            vwind300_data = extract_data(vwind_file, vwind_variable_name, 300, time)
        elif(949 >= tc_central_pressure[count] >= 940):
            uwind850_data = extract_data(uwind_file, uwind_variable_name, 850, time)
            vwind850_data = extract_data(vwind_file, vwind_variable_name, 850, time)
            uwind800_data = extract_data(uwind_file, uwind_variable_name, 800, time)
            vwind800_data = extract_data(vwind_file, vwind_variable_name, 800, time)
            uwind750_data = extract_data(uwind_file, uwind_variable_name, 750, time)
            vwind750_data = extract_data(vwind_file, vwind_variable_name, 750, time)
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
            uwind650_data = extract_data(uwind_file, uwind_variable_name, 650, time)
            vwind650_data = extract_data(vwind_file, vwind_variable_name, 650, time)
            uwind600_data = extract_data(uwind_file, uwind_variable_name, 600, time)
            vwind600_data = extract_data(vwind_file, vwind_variable_name, 600, time)
            uwind550_data = extract_data(uwind_file, uwind_variable_name, 550, time)
            vwind550_data = extract_data(vwind_file, vwind_variable_name, 550, time)
            uwind500_data = extract_data(uwind_file, uwind_variable_name, 500, time)
            vwind500_data = extract_data(vwind_file, vwind_variable_name, 500, time)
            uwind450_data = extract_data(uwind_file, uwind_variable_name, 450, time)
            vwind450_data = extract_data(vwind_file, vwind_variable_name, 450, time)
            uwind400_data = extract_data(uwind_file, uwind_variable_name, 400, time)
            vwind400_data = extract_data(vwind_file, vwind_variable_name, 400, time)
            uwind350_data = extract_data(uwind_file, uwind_variable_name, 350, time)
            vwind350_data = extract_data(vwind_file, vwind_variable_name, 350, time)
            uwind300_data = extract_data(uwind_file, uwind_variable_name, 300, time)
            vwind300_data = extract_data(vwind_file, vwind_variable_name, 300, time)
            uwind250_data = extract_data(uwind_file, uwind_variable_name, 250, time)
            vwind250_data = extract_data(vwind_file, vwind_variable_name, 250, time)
        else:
            uwind700_data = extract_data(uwind_file, uwind_variable_name, 700, time)
            vwind700_data = extract_data(vwind_file, vwind_variable_name, 700, time)
            uwind650_data = extract_data(uwind_file, uwind_variable_name, 650, time)
            vwind650_data = extract_data(vwind_file, vwind_variable_name, 650, time)
            uwind600_data = extract_data(uwind_file, uwind_variable_name, 600, time)
            vwind600_data = extract_data(vwind_file, vwind_variable_name, 600, time)
            uwind550_data = extract_data(uwind_file, uwind_variable_name, 550, time)
            vwind550_data = extract_data(vwind_file, vwind_variable_name, 550, time)
            uwind500_data = extract_data(uwind_file, uwind_variable_name, 500, time)
            vwind500_data = extract_data(vwind_file, vwind_variable_name, 500, time)
            uwind450_data = extract_data(uwind_file, uwind_variable_name, 450, time)
            vwind450_data = extract_data(vwind_file, vwind_variable_name, 450, time)
            uwind400_data = extract_data(uwind_file, uwind_variable_name, 400, time)
            vwind400_data = extract_data(vwind_file, vwind_variable_name, 400, time)
            uwind350_data = extract_data(uwind_file, uwind_variable_name, 350, time)
            vwind350_data = extract_data(vwind_file, vwind_variable_name, 350, time)
            uwind300_data = extract_data(uwind_file, uwind_variable_name, 300, time)
            vwind300_data = extract_data(vwind_file, vwind_variable_name, 300, time)
            uwind250_data = extract_data(uwind_file, uwind_variable_name, 250, time)
            vwind250_data = extract_data(vwind_file, vwind_variable_name, 250, time)
            uwind200_data = extract_data(uwind_file, uwind_variable_name, 200, time)
            vwind200_data = extract_data(vwind_file, vwind_variable_name, 200, time)
        
        # Defining the steering layer: Dividing the steering layer into sublayers with
        # dp = 50 hPa thickness, i.e. into pairs of pressure levels, and concatenating
        # the wind data in each of the included pairs of pressure levels.
        # NOTE: The extracted pressure levels from ERA5 dataset are in
        # constant increments of dp = 50 hPa.
        if(tc_central_pressure[count] >= 1000):
            steering_layer = 'Steering layer: 850–700 hPa'
            DP = float(850-700) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
        elif(999 >= tc_central_pressure[count] >= 990):
            steering_layer = 'Steering layer: 850–500 hPa'
            DP = float(850-500) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            ulayer_700to650 = xr.concat([uwind700_data, uwind650_data], dim='pressure_level')
            ulayer_650to600 = xr.concat([uwind650_data, uwind600_data], dim='pressure_level')
            ulayer_600to550 = xr.concat([uwind600_data, uwind550_data], dim='pressure_level')
            ulayer_550to500 = xr.concat([uwind550_data, uwind500_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
            vlayer_700to650 = xr.concat([vwind700_data, vwind650_data], dim='pressure_level')
            vlayer_650to600 = xr.concat([vwind650_data, vwind600_data], dim='pressure_level')
            vlayer_600to550 = xr.concat([vwind600_data, vwind550_data], dim='pressure_level')
            vlayer_550to500 = xr.concat([vwind550_data, vwind500_data], dim='pressure_level')
        elif(989 >= tc_central_pressure[count] >= 970):
            steering_layer = 'Steering layer: 850–400 hPa'
            DP = float(850-400) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            ulayer_700to650 = xr.concat([uwind700_data, uwind650_data], dim='pressure_level')
            ulayer_650to600 = xr.concat([uwind650_data, uwind600_data], dim='pressure_level')
            ulayer_600to550 = xr.concat([uwind600_data, uwind550_data], dim='pressure_level')
            ulayer_550to500 = xr.concat([uwind550_data, uwind500_data], dim='pressure_level')
            ulayer_500to450 = xr.concat([uwind500_data, uwind450_data], dim='pressure_level')
            ulayer_450to400 = xr.concat([uwind450_data, uwind400_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
            vlayer_700to650 = xr.concat([vwind700_data, vwind650_data], dim='pressure_level')
            vlayer_650to600 = xr.concat([vwind650_data, vwind600_data], dim='pressure_level')
            vlayer_600to550 = xr.concat([vwind600_data, vwind550_data], dim='pressure_level')
            vlayer_550to500 = xr.concat([vwind550_data, vwind500_data], dim='pressure_level')
            vlayer_500to450 = xr.concat([vwind500_data, vwind450_data], dim='pressure_level')
            vlayer_450to400 = xr.concat([vwind450_data, vwind400_data], dim='pressure_level')
        elif(969 >= tc_central_pressure[count] >= 950):
            steering_layer = 'Steering layer: 850–300 hPa'
            DP = float(850-300) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            ulayer_700to650 = xr.concat([uwind700_data, uwind650_data], dim='pressure_level')
            ulayer_650to600 = xr.concat([uwind650_data, uwind600_data], dim='pressure_level')
            ulayer_600to550 = xr.concat([uwind600_data, uwind550_data], dim='pressure_level')
            ulayer_550to500 = xr.concat([uwind550_data, uwind500_data], dim='pressure_level')
            ulayer_500to450 = xr.concat([uwind500_data, uwind450_data], dim='pressure_level')
            ulayer_450to400 = xr.concat([uwind450_data, uwind400_data], dim='pressure_level')
            ulayer_400to350 = xr.concat([uwind400_data, uwind350_data], dim='pressure_level')
            ulayer_350to300 = xr.concat([uwind350_data, uwind300_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
            vlayer_700to650 = xr.concat([vwind700_data, vwind650_data], dim='pressure_level')
            vlayer_650to600 = xr.concat([vwind650_data, vwind600_data], dim='pressure_level')
            vlayer_600to550 = xr.concat([vwind600_data, vwind550_data], dim='pressure_level')
            vlayer_550to500 = xr.concat([vwind550_data, vwind500_data], dim='pressure_level')
            vlayer_500to450 = xr.concat([vwind500_data, vwind450_data], dim='pressure_level')
            vlayer_450to400 = xr.concat([vwind450_data, vwind400_data], dim='pressure_level')
            vlayer_400to350 = xr.concat([vwind400_data, vwind350_data], dim='pressure_level')
            vlayer_350to300 = xr.concat([vwind350_data, vwind300_data], dim='pressure_level')
        elif(949 >= tc_central_pressure[count] >= 940):
            steering_layer = 'Steering layer: 850–250 hPa'
            DP = float(850-250) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            ulayer_700to650 = xr.concat([uwind700_data, uwind650_data], dim='pressure_level')
            ulayer_650to600 = xr.concat([uwind650_data, uwind600_data], dim='pressure_level')
            ulayer_600to550 = xr.concat([uwind600_data, uwind550_data], dim='pressure_level')
            ulayer_550to500 = xr.concat([uwind550_data, uwind500_data], dim='pressure_level')
            ulayer_500to450 = xr.concat([uwind500_data, uwind450_data], dim='pressure_level')
            ulayer_450to400 = xr.concat([uwind450_data, uwind400_data], dim='pressure_level')
            ulayer_400to350 = xr.concat([uwind400_data, uwind350_data], dim='pressure_level')
            ulayer_350to300 = xr.concat([uwind350_data, uwind300_data], dim='pressure_level')
            ulayer_300to250 = xr.concat([uwind300_data, uwind250_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
            vlayer_700to650 = xr.concat([vwind700_data, vwind650_data], dim='pressure_level')
            vlayer_650to600 = xr.concat([vwind650_data, vwind600_data], dim='pressure_level')
            vlayer_600to550 = xr.concat([vwind600_data, vwind550_data], dim='pressure_level')
            vlayer_550to500 = xr.concat([vwind550_data, vwind500_data], dim='pressure_level')
            vlayer_500to450 = xr.concat([vwind500_data, vwind450_data], dim='pressure_level')
            vlayer_450to400 = xr.concat([vwind450_data, vwind400_data], dim='pressure_level')
            vlayer_400to350 = xr.concat([vwind400_data, vwind350_data], dim='pressure_level')
            vlayer_350to300 = xr.concat([vwind350_data, vwind300_data], dim='pressure_level')
            vlayer_300to250 = xr.concat([vwind300_data, vwind250_data], dim='pressure_level')
        else:
            steering_layer = 'Steering layer: 700–200 hPa'
            DP = float(700-200) # Total layer thickness in hPa (difference of vertical pressure values)
            ulayer_850to800 = xr.concat([uwind850_data, uwind800_data], dim='pressure_level')
            ulayer_800to750 = xr.concat([uwind800_data, uwind750_data], dim='pressure_level')
            ulayer_750to700 = xr.concat([uwind750_data, uwind700_data], dim='pressure_level')
            ulayer_700to650 = xr.concat([uwind700_data, uwind650_data], dim='pressure_level')
            ulayer_650to600 = xr.concat([uwind650_data, uwind600_data], dim='pressure_level')
            ulayer_600to550 = xr.concat([uwind600_data, uwind550_data], dim='pressure_level')
            ulayer_550to500 = xr.concat([uwind550_data, uwind500_data], dim='pressure_level')
            ulayer_500to450 = xr.concat([uwind500_data, uwind450_data], dim='pressure_level')
            ulayer_450to400 = xr.concat([uwind450_data, uwind400_data], dim='pressure_level')
            ulayer_400to350 = xr.concat([uwind400_data, uwind350_data], dim='pressure_level')
            ulayer_350to300 = xr.concat([uwind350_data, uwind300_data], dim='pressure_level')
            ulayer_300to250 = xr.concat([uwind300_data, uwind250_data], dim='pressure_level')
            ulayer_250to200 = xr.concat([uwind250_data, uwind200_data], dim='pressure_level')
            vlayer_850to800 = xr.concat([vwind850_data, vwind800_data], dim='pressure_level')
            vlayer_800to750 = xr.concat([vwind800_data, vwind750_data], dim='pressure_level')
            vlayer_750to700 = xr.concat([vwind750_data, vwind700_data], dim='pressure_level')
            vlayer_700to650 = xr.concat([vwind700_data, vwind650_data], dim='pressure_level')
            vlayer_650to600 = xr.concat([vwind650_data, vwind600_data], dim='pressure_level')
            vlayer_600to550 = xr.concat([vwind600_data, vwind550_data], dim='pressure_level')
            vlayer_550to500 = xr.concat([vwind550_data, vwind500_data], dim='pressure_level')
            vlayer_500to450 = xr.concat([vwind500_data, vwind450_data], dim='pressure_level')
            vlayer_450to400 = xr.concat([vwind450_data, vwind400_data], dim='pressure_level')
            vlayer_400to350 = xr.concat([vwind400_data, vwind350_data], dim='pressure_level')
            vlayer_350to300 = xr.concat([vwind350_data, vwind300_data], dim='pressure_level')
            vlayer_300to250 = xr.concat([vwind300_data, vwind250_data], dim='pressure_level')
            vlayer_250to200 = xr.concat([vwind250_data, vwind200_data], dim='pressure_level')
        
        # Computing the components of the deep-layer mean wind
        if(tc_central_pressure[count] >= 1000):
            u_steeringwind = ((ulayer_850to800.mean(dim='pressure_level')*dp) + (ulayer_800to750.mean(dim='pressure_level')*dp) + (ulayer_750to700.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_850to800.mean(dim='pressure_level')*dp) + (vlayer_800to750.mean(dim='pressure_level')*dp) + (vlayer_750to700.mean(dim='pressure_level')*dp)) / DP
        elif(999 >= tc_central_pressure[count] >= 990):
            u_steeringwind = ((ulayer_850to800.mean(dim='pressure_level')*dp) + (ulayer_800to750.mean(dim='pressure_level')*dp) + (ulayer_750to700.mean(dim='pressure_level')*dp) + (ulayer_700to650.mean(dim='pressure_level')*dp) + (ulayer_650to600.mean(dim='pressure_level')*dp) + (ulayer_600to550.mean(dim='pressure_level')*dp) + (ulayer_550to500.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_850to800.mean(dim='pressure_level')*dp) + (vlayer_800to750.mean(dim='pressure_level')*dp) + (vlayer_750to700.mean(dim='pressure_level')*dp) + (vlayer_700to650.mean(dim='pressure_level')*dp) + (vlayer_650to600.mean(dim='pressure_level')*dp) + (vlayer_600to550.mean(dim='pressure_level')*dp) + (vlayer_550to500.mean(dim='pressure_level')*dp)) / DP
        elif(989 >= tc_central_pressure[count] >= 970):
            u_steeringwind = ((ulayer_850to800.mean(dim='pressure_level')*dp) + (ulayer_800to750.mean(dim='pressure_level')*dp) + (ulayer_750to700.mean(dim='pressure_level')*dp) + (ulayer_700to650.mean(dim='pressure_level')*dp) + (ulayer_650to600.mean(dim='pressure_level')*dp) + (ulayer_600to550.mean(dim='pressure_level')*dp) + (ulayer_550to500.mean(dim='pressure_level')*dp) + (ulayer_500to450.mean(dim='pressure_level')*dp) + (ulayer_450to400.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_850to800.mean(dim='pressure_level')*dp) + (vlayer_800to750.mean(dim='pressure_level')*dp) + (vlayer_750to700.mean(dim='pressure_level')*dp) + (vlayer_700to650.mean(dim='pressure_level')*dp) + (vlayer_650to600.mean(dim='pressure_level')*dp) + (vlayer_600to550.mean(dim='pressure_level')*dp) + (vlayer_550to500.mean(dim='pressure_level')*dp) + (vlayer_500to450.mean(dim='pressure_level')*dp) + (vlayer_450to400.mean(dim='pressure_level')*dp)) / DP
        elif(969 >= tc_central_pressure[count] >= 950):
            u_steeringwind = ((ulayer_850to800.mean(dim='pressure_level')*dp) + (ulayer_800to750.mean(dim='pressure_level')*dp) + (ulayer_750to700.mean(dim='pressure_level')*dp) + (ulayer_700to650.mean(dim='pressure_level')*dp) + (ulayer_650to600.mean(dim='pressure_level')*dp) + (ulayer_600to550.mean(dim='pressure_level')*dp) + (ulayer_550to500.mean(dim='pressure_level')*dp) + (ulayer_500to450.mean(dim='pressure_level')*dp) + (ulayer_450to400.mean(dim='pressure_level')*dp) + (ulayer_400to350.mean(dim='pressure_level')*dp) + (ulayer_350to300.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_850to800.mean(dim='pressure_level')*dp) + (vlayer_800to750.mean(dim='pressure_level')*dp) + (vlayer_750to700.mean(dim='pressure_level')*dp) + (vlayer_700to650.mean(dim='pressure_level')*dp) + (vlayer_650to600.mean(dim='pressure_level')*dp) + (vlayer_600to550.mean(dim='pressure_level')*dp) + (vlayer_550to500.mean(dim='pressure_level')*dp) + (vlayer_500to450.mean(dim='pressure_level')*dp) + (vlayer_450to400.mean(dim='pressure_level')*dp) + (vlayer_400to350.mean(dim='pressure_level')*dp) + (vlayer_350to300.mean(dim='pressure_level')*dp)) / DP
        elif(949 >= tc_central_pressure[count] >= 940):
            u_steeringwind = ((ulayer_850to800.mean(dim='pressure_level')*dp) + (ulayer_800to750.mean(dim='pressure_level')*dp) + (ulayer_750to700.mean(dim='pressure_level')*dp) + (ulayer_700to650.mean(dim='pressure_level')*dp) + (ulayer_650to600.mean(dim='pressure_level')*dp) + (ulayer_600to550.mean(dim='pressure_level')*dp) + (ulayer_550to500.mean(dim='pressure_level')*dp) + (ulayer_500to450.mean(dim='pressure_level')*dp) + (ulayer_450to400.mean(dim='pressure_level')*dp) + (ulayer_400to350.mean(dim='pressure_level')*dp) + (ulayer_350to300.mean(dim='pressure_level')*dp) + (ulayer_300to250.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_850to800.mean(dim='pressure_level')*dp) + (vlayer_800to750.mean(dim='pressure_level')*dp) + (vlayer_750to700.mean(dim='pressure_level')*dp) + (vlayer_700to650.mean(dim='pressure_level')*dp) + (vlayer_650to600.mean(dim='pressure_level')*dp) + (vlayer_600to550.mean(dim='pressure_level')*dp) + (vlayer_550to500.mean(dim='pressure_level')*dp) + (vlayer_500to450.mean(dim='pressure_level')*dp) + (vlayer_450to400.mean(dim='pressure_level')*dp) + (vlayer_400to350.mean(dim='pressure_level')*dp) + (vlayer_350to300.mean(dim='pressure_level')*dp) + (vlayer_300to250.mean(dim='pressure_level')*dp)) / DP
        else:
            u_steeringwind = ((ulayer_700to650.mean(dim='pressure_level')*dp) + (ulayer_650to600.mean(dim='pressure_level')*dp) + (ulayer_600to550.mean(dim='pressure_level')*dp) + (ulayer_550to500.mean(dim='pressure_level')*dp) + (ulayer_500to450.mean(dim='pressure_level')*dp) + (ulayer_450to400.mean(dim='pressure_level')*dp) + (ulayer_400to350.mean(dim='pressure_level')*dp) + (ulayer_350to300.mean(dim='pressure_level')*dp) + (ulayer_300to250.mean(dim='pressure_level')*dp) + (ulayer_250to200.mean(dim='pressure_level')*dp)) / DP
            v_steeringwind = ((vlayer_700to650.mean(dim='pressure_level')*dp) + (vlayer_650to600.mean(dim='pressure_level')*dp) + (vlayer_600to550.mean(dim='pressure_level')*dp) + (vlayer_550to500.mean(dim='pressure_level')*dp) + (vlayer_500to450.mean(dim='pressure_level')*dp) + (vlayer_450to400.mean(dim='pressure_level')*dp) + (vlayer_400to350.mean(dim='pressure_level')*dp) + (vlayer_350to300.mean(dim='pressure_level')*dp) + (vlayer_300to250.mean(dim='pressure_level')*dp) + (vlayer_250to200.mean(dim='pressure_level')*dp)) / DP
        
        # Computing speed of the deep-layer mean wind (magnitude of the wind velocity vectors) and relative vorticity
        if(product == 'Wind velocity'):
            steering_wind_speed = steering_wind_speed_calc(u_steeringwind, v_steeringwind) # In km/h
        else:
            rel_vorticity = mpcalc.vorticity(u_steeringwind*units('m/s'), v_steeringwind*units('m/s')) # In 1/s
        
        # Initializing the map
        if(product == 'Wind velocity'):
            x, y = np.meshgrid(steering_wind_speed['longitude'], steering_wind_speed['latitude'])
        else:
            x, y = np.meshgrid(rel_vorticity['longitude'], rel_vorticity['latitude'])
            
        fig = plt.figure(figsize=(25., 25.), dpi=250)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([90, 180, -10, 50], ccrs.PlateCarree())
        # ax.set_extent([100, 150, -5, 40], ccrs.PlateCarree())
        
        # Add gridlines
        gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 10), ylocs=np.arange(-90, 91, 10), color='gray', linestyle='--')
        # gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 5), ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
        gridlines.top_labels = False
        gridlines.right_labels = False
        
        # Adding the boundary of the Philippine Area of Responsibility (PAR)
        PAR_longitudes = [120, 135, 135, 115, 115, 120, 120]
        PAR_latitudes = [25, 25, 5, 5, 15, 21, 25]
        ax.plot(PAR_longitudes, PAR_latitudes, color='black', linewidth=3, transform=ccrs.PlateCarree(), linestyle='--')
        
        # Visualizing wind speed and magnitude of relative vorticity
        if(product == 'Wind velocity'):
            wind_speed_visuals = ax.contourf(x, y, steering_wind_speed, transform=ccrs.PlateCarree(), levels=np.arange(0, 200, 1), cmap=windspeed_colormap)
        else:
            rel_vorticity_visuals = ax.contourf(x, y, rel_vorticity, transform=ccrs.PlateCarree(), levels=np.arange(-0.0009, 0.0009, 0.00001), cmap=vorticity_colormap)
        
        # Adding the coastlines
        ax.coastlines('10m', edgecolor='black', linewidth=2.5)
        
        # Adding streamlines for wind direction
        streamlines = ax.streamplot(x, y, u_steeringwind, v_steeringwind, color='k', density=2.5, arrowsize=3, arrowstyle='->')
        
        # Plotting the TC's current location and 6-hour displacement
        latitude_fordisplacement = [latitude[count], latitude[count+1]]
        longitude_fordisplacement = [longitude[count], longitude[count+1]]
        latitude_forcurrentposition = latitude[count]
        longitude_forcurrentposition = longitude[count]
        
        ax.plot(longitude_fordisplacement, latitude_fordisplacement, 'k-X', transform=ccrs.PlateCarree(), linewidth=5, markersize=15)
        ax.plot(longitude_forcurrentposition, latitude_forcurrentposition, 'ro', transform=ccrs.PlateCarree(), markersize=15)
        
        # Calculating the TC's average forward speed
        current_position = (latitude[count], longitude[count])
        position_6hrlater = (latitude[count+1], longitude[count+1])
        
        displacement_6hr = geodesic(current_position, position_6hrlater)
        forward_speed = round(displacement_6hr.kilometers / 6, 1)
        
        # Adding a colorbar
        if(product == 'Wind velocity'):
            cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
            norm = Normalize(vmin=0, vmax=200)
            colorbar = ColorbarBase(cax, cmap=windspeed_colormap, norm=norm, extend='max', ticks=np.arange(0, 201, 5))
            colorbar.set_label('Wind speed (km/h)')
        else:
            cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
            norm = Normalize(vmin=-0.0009, vmax=0.0009)
            colorbar = ColorbarBase(cax, cmap=vorticity_colormap, norm=norm, extend='both', ticks=np.arange(-0.0009, 0.00091, 0.00009))
            colorbar.set_label('Relative vorticity (1/s)')
        
        # Adding plot title and captions
        title = 'Environmental steering [' + time[:10] + ' ' + time[11:16] + ' UTC, ' + steering_layer + ']\nTyphoon SAOLA @ CentPress: ' + str(tc_central_pressure[count]) + ' hPa, 10-min MSW: ' + str(round((max_sustained_winds_knots[count]*1.852)/5)*5) + ' km/h, Ave. forward speed: ' + str(forward_speed) + ' km/h'
        plt.suptitle(title, x=0.125, y=0.786, fontsize=20, ha='left', va='top')
        # plt.suptitle(title, x=0.125, y=0.877, fontsize=20, ha='left', va='top')
        
        plt.show()
    
        # Printing the status of rendering the visualization
        if(product == 'Wind velocity'):
            print('\n', product,'>> Slide', count+1, 'DONE!')
        else:
            print('\n', product,'>> Slide', count+1, 'DONE!')
        
        count = count + 1
        
        # Should stop at index = 37, this is to make sure the program ends without incurring any errors
        if(count == 37):
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
