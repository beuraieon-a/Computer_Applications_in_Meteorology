# Program for visualizing visible satellite imagery from Himawari-8/9 as a
# "true color" composite of the veggie near-IR band and the blue and red visible bands,
# with a nighttime view using an RGB composite of near-IR (shortwave) and
# clean IR longwave bands

import os
import datetime
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Satellite data extraction for a specific area
def extract_satdata(file_path, band, lat_bounds, lon_bounds):
    extracted_satdata = xr.open_dataset(file_path)[band].sel(latitude=lat_bounds, longitude=lon_bounds)
    return extracted_satdata

# FUNCTION 2: Enhancement of color saturation
def enhance_saturation(image, saturation_factor):
    mean = np.mean(image, axis=-1, keepdims=True)
    return np.clip(mean + (image - mean) * saturation_factor, 0, 1)

# FUNCTION 3: Application of gamma correction
def apply_gamma_correction(image, gamma):
    return np.clip(image ** gamma, 0, 1)

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Setting the data extraction parameters
lat_south = 0
lat_north = 30
lon_west = 108
lon_east = 146

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
    
    #### Extracting satellite imagery data ####

    # For this visualization, the following bands are needed:
    # Band 3 - Visible, red (0.64 µm)
    # Band 4 - Near-infrared, veggie (0.86 µm), as a substitute to the original green visible band to better represent landmasses
    # Band 1: Visible, blue (0.47 µm)
    # Band 7: Infrared, shortwave (3.9 µm)
    # Band 13: Infrared, clean IR longwave (10.4 µm)
    # The program will also extract data on solar zenith angle (the angle between the sun and the local vertical / zenith).

    VIS_red = extract_satdata(sat_path, 'albedo_03', lat_boundary, lon_boundary)
    NIR_green = extract_satdata(sat_path, 'albedo_04', lat_boundary, lon_boundary)
    VIS_blue = extract_satdata(sat_path, 'albedo_01', lat_boundary, lon_boundary)
    IR_shortwave = extract_satdata(sat_path, 'tbb_07', lat_boundary, lon_boundary)
    IR_cloudtop = extract_satdata(sat_path, 'tbb_13', lat_boundary, lon_boundary)  
    solar_zenith_angle = extract_satdata(sat_path, 'SOZ', lat_boundary, lon_boundary)

    # Extracting the observation time
    time = str(xr.open_dataset(sat_path)['start_time'].values[0])

    #### Preparing the daytime imagery ####

    # Creating the visible+near-IR RGB channels by normalizing the data to the range 0-1
    VIS_red = np.clip(VIS_red, 0, 1)
    NIR_green = np.clip(NIR_green, 0, 1)
    VIS_blue = np.clip(VIS_blue, 0, 1)

    # Toning down the NIR_green band to make it a "true" green (making the image not too vibrant)
    # by interpolating the value to simulate a natural green color
    true_green = (0.45 * VIS_red) + (0.1 * NIR_green) + (0.45 * VIS_blue)
    true_green = np.clip(true_green, 0, 1)  # Apply 0-1 normalization, just in case

    # Creating the "true color" RGB composite
    true_color = np.dstack([VIS_red, true_green, VIS_blue])
    true_color = enhance_saturation(true_color, 2.0)
    true_color = apply_gamma_correction(true_color, 0.5)

    #### Preparing the nighttime imagery ####

    # Normalizing the IR datasets within a range of brightness temperatures in kelvin
    # Note: 183.15 K = -90°C, 193.15 K = -80°C, 313.15 K = 40°C
    IR_shortwave = (IR_shortwave - 193.15) / (313.15 - 193.15)
    IR_cloudtop = (IR_cloudtop - 183.15) / (313.15 - 183.15)

    # Creating the IR RGB channels by normalizing the data to the range 0-1
    night_red = np.clip(IR_shortwave, 0, 1) # Red: Band 7 (Infrared, shortwave)
    night_green = np.clip(IR_shortwave, 0, 1) # Green: Band 7 (Infrared, shortwave)
    night_blue = np.clip(IR_cloudtop, 0, 1) # Blue: Band 13 (Infrared, clean IR longwave)

    # Inverting the IR colors so that cold cloud-tops are white
    night_red = 1 - night_red
    night_green = 1 - night_green
    night_blue = 1 - night_blue

    # Increasing the brightness of the IR shortwave band
    night_red = apply_gamma_correction(night_red, 0.98)
    night_green = apply_gamma_correction(night_green, 0.98)

    # Creating the night-side RGB composite imagery
    night_view = np.dstack([night_red, night_green, night_blue])
    night_view = enhance_saturation(night_view, 2.0)

    #### Preparing the day-night composite imagery ####

    # Creating a day-night mask based on the solar zenith angle
    day_mask = np.ones_like(solar_zenith_angle)

    # >> Gradually transitions to night view for grid points near the terminator,
    # i.e. solar zenith angle = 80° to 90°
    day_mask = np.where(solar_zenith_angle >= 80,
                        np.clip(1 - (solar_zenith_angle - 80) / 10, 0, 1),
                        day_mask) 

    # Day view completely transparent, i.e. night view only,
    # for grid points with solar zenith angle > 90°
    day_mask = np.where(solar_zenith_angle > 90, 0, day_mask)  

    # Expanding the masks to match the RGB dimensions
    day_mask_rgb = np.dstack([day_mask] * 3)
    night_mask_rgb = 1 - day_mask_rgb

    # Blending the day and night images using the masks
    daynight_composite = (true_color * day_mask_rgb) + (night_view * night_mask_rgb)

    #### Visualizing the satellite imagery ####

    # Initializing the map
    x, y = np.meshgrid(VIS_red['longitude'], VIS_red['latitude'])
    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], ccrs.PlateCarree())
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Plotting the satellite imagery
    satellite_image = ax.imshow(daynight_composite, origin='upper',
                                extent=[lon_west, lon_east, lat_south, lat_north],
                                transform=ccrs.PlateCarree())

    # Adding the coastlines
    ax.coastlines('10m', edgecolor='black', linewidth=2)

    # Adding gridlines
    gridlines = ax.gridlines(draw_labels=False, xlocs=np.arange(-180, 181, 2),
                             ylocs=np.arange(-90, 91, 2), color='gray',
                             linestyle='--')
    xticks = np.arange(lon_west, lon_east + 1, 2)
    yticks = np.arange(lat_south, lat_north + 1, 2)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    # Plotting the TCID boundary
    TCID, = ax.plot(TCID_longitudes, TCID_latitudes, color='yellow', linewidth=3,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Plotting the TCAD boundary
    TCAD, = ax.plot(TCAD_longitudes, TCAD_latitudes, color='orange', linewidth=3,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Plotting the PAR boundary
    PAR, = ax.plot(PAR_longitudes, PAR_latitudes, color='red', linewidth=3,
            transform=ccrs.PlateCarree(), linestyle='-')

    # Adding plot title and captions
    time_truncated = time[:16]
    time_obj = datetime.datetime.strptime(time_truncated, '%Y-%m-%dT%H:%M')
    formatted_date = time_obj.strftime('%d %B %Y')
    formatted_time = time_obj.strftime('%H:%M')

    title = 'Satellite image: True color RGB (day) and infrared (night) | Himawari-9'
    caption = f'{formatted_time} UTC, {formatted_date}'
    plt.suptitle(title, x=0.125, y=0.844, fontsize=35, ha='left', va='top')
    plt.figtext(0.125, 0.821, caption, fontsize=29, ha='left', va='top')

    # Showing the visualized satellite imagery
    plt.show()

# End of program