# Program for determining the number of hours
# before a specific grid location attains a given wind speed threshold

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction for wind data
def extract_data(nc_file, variable_name, start_time, end_time):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(time=slice(start_time, end_time))
    dataset.close()

    return extracted_data

# FUNCTION 2: Data extraction for time data
def extract_time(nc_file, variable_name, start_time, end_time):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(time=slice(start_time, end_time)).values
    extracted_data_list = extracted_data.tolist()
    
    extracted_data_list = [np.datetime64(int(dt), 'ns') for dt in extracted_data_list]
    datetime_str_list = [str(dt) for dt in extracted_data_list]
    
    dataset.close()

    return datetime_str_list

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

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Opening files
uwind_file = 'C:/Users/Brian/Desktop/Future thesis works/VAR_100U.nc'
vwind_file = 'C:/Users/Brian/Desktop/Future thesis works/VAR_100V.nc'
TCWS_rgb_file_path = 'C:/Users/Brian/Desktop/Future thesis works/tcws.rgb'

# Colormap
TCWS_colors = read_rgb_file(TCWS_rgb_file_path)
TCWS_colormap = mcolors.ListedColormap(TCWS_colors)

# Preparing wind variable names as inputs to FUNCTION 1 (extract_data)
uwind_variable_name = 'VAR_100U'
vwind_variable_name = 'VAR_100V'
time_variable_name = 'time'

# Data extraction
start_time = '2022-09-24T05:00:00'
end_time = '2022-09-26T05:00:00'

time_data = extract_time(uwind_file, time_variable_name, start_time, end_time)
uwind_data = extract_data(uwind_file, uwind_variable_name, start_time, end_time)
vwind_data = extract_data(vwind_file, vwind_variable_name, start_time, end_time)

# Wind speed calculation
wind_speed = steering_wind_speed_calc(uwind_data, vwind_data)

lat_min, lat_max = 4, 21.5
lon_min, lon_max = 114, 127

TD_wind_threshold = 39
TS_wind_threshold = 62

TD_time_elapsed = xr.DataArray(np.full((len(uwind_data.latitude), len(uwind_data.longitude)), np.nan),
                            coords=[uwind_data.latitude, uwind_data.longitude],
                            dims=['latitude', 'longitude'])

TS_time_elapsed = xr.DataArray(np.full((len(uwind_data.latitude), len(uwind_data.longitude)), np.nan),
                            coords=[uwind_data.latitude, uwind_data.longitude],
                            dims=['latitude', 'longitude'])

TCWS = xr.DataArray(np.full((len(uwind_data.latitude), len(uwind_data.longitude)), np.nan),
                            coords=[uwind_data.latitude, uwind_data.longitude],
                            dims=['latitude', 'longitude'])

print('\nDetermining time elapsed before arrival of TD winds...\n')

for lat in np.arange(lat_min, lat_max, 0.25):
    for lon in np.arange(lon_min, lon_max, 0.25):
        count = 0
        for time in time_data:
            wind_speed_at_latlon = wind_speed.sel(latitude=lat, longitude=lon, time=time, method='nearest').item()

            if wind_speed_at_latlon >= TD_wind_threshold:
                break
            else:
                count = count + 1

        # Update the time-elapsed value for the current latitude and longitude
        TD_time_elapsed.loc[dict(latitude=lat, longitude=lon)] = count
        print('\nDone analyzing at:', lat, lon)

print('\nDetermining time elapsed before arrival of TS winds...\n')

for lat in np.arange(lat_min, lat_max, 0.25):
    for lon in np.arange(lon_min, lon_max, 0.25):
        count = 0
        for time in time_data:
            wind_speed_at_latlon = wind_speed.sel(latitude=lat, longitude=lon, time=time, method='nearest').item()

            if wind_speed_at_latlon >= TS_wind_threshold:
                break
            else:
                count = count + 1

        # Update the time-elapsed value for the current latitude and longitude
        TS_time_elapsed.loc[dict(latitude=lat, longitude=lon)] = count
        print('\nDone analyzing at:', lat, lon)

print('\nDetermining TCWS...\n')

for lat in np.arange(lat_min, lat_max, 0.25):
    for lon in np.arange(lon_min, lon_max, 0.25):
        time_elapsed_TD_at_latlon  = TD_time_elapsed.sel(latitude=lat, longitude=lon, method='nearest').item()
        time_elapsed_TS_at_latlon  = TS_time_elapsed.sel(latitude=lat, longitude=lon, method='nearest').item()
        
        if time_elapsed_TS_at_latlon <= 24:
            TCWS.loc[dict(latitude=lat, longitude=lon)] = 2
        elif time_elapsed_TD_at_latlon <= 36:
            TCWS.loc[dict(latitude=lat, longitude=lon)] = 1
        else:
            TCWS.loc[dict(latitude=lat, longitude=lon)] = 0

        print('\nDone analyzing at:', lat, lon)
        print('time_elapsed_TD_at_latlon:', time_elapsed_TD_at_latlon)
        print('time_elapsed_TS_at_latlon:', time_elapsed_TS_at_latlon)
        print('TCWS:', TCWS.loc[dict(latitude=lat, longitude=lon)].item())

# Visualization

print(TCWS)

x, y = np.meshgrid(uwind_data['longitude'], uwind_data['latitude'])
fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([114, 130, 4, 21.5], ccrs.PlateCarree())

gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 1), ylocs=np.arange(-90, 91, 1), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

TCWS_visuals = ax.contourf(x, y, TCWS, transform=ccrs.PlateCarree(), levels=np.arange(0, 5.01, 0.01), cmap=TCWS_colormap)
ax.coastlines('10m', edgecolor='black', linewidth=2.5)

TCWS.to_netcdf('C:/Users/Brian/Desktop/TCWS.nc')

# print(TD_time_elapsed)

# End of program
