# Program for extracting TC track data from a plain text (TXT) file

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Opening file
tctrack_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/saola_track.txt'

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

tctrack_data_dict = {
    'Date & Time (UTC)': time_frame,
    'TC Category': tc_category,
    'Latitude': latitude,
    'Longitude': longitude,
    'TC Central Pressure (hPa)': tc_central_pressure,
    'Max Sustained Winds (kts)': max_sustained_winds_knots
}

tctrack_data_df = pd.DataFrame(tctrack_data_dict)

fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([116, 130, 14, 24], ccrs.PlateCarree())

gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 2), ylocs=np.arange(-90, 91, 2), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

ax.coastlines('10m', edgecolor='black', linewidth=2.5)
ax.plot(tctrack_data_df['Longitude'], tctrack_data_df['Latitude'], '-k', transform=ccrs.PlateCarree(), markersize=15)
ax.plot(tctrack_data_df['Longitude'], tctrack_data_df['Latitude'], 'or', transform=ccrs.PlateCarree(), markersize=15)

plt.show()