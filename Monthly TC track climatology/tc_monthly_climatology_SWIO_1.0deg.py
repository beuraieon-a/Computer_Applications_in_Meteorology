# Program for generating monthly climatology maps of
# Southwest Indian tropical cyclone tracks using IBTrACS data

import numpy as np
import pandas as pd
from scipy.stats import mode
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Extracting TC best track data from the IBTrACS CSV files
all_TC_tracks = pd.read_csv('C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/ibtracs.SI.list.v04r00.csv',
                      usecols=['SID', 'ISO_TIME', 'LAT', 'LON'])

# Converting the string times into datetime type
all_TC_tracks['ISO_TIME'] = pd.to_datetime(all_TC_tracks['ISO_TIME'])

# Extracting the month value from the datetimes
all_TC_tracks['MONTH'] = all_TC_tracks['ISO_TIME'].dt.month

# Determining the month of occurence of each TC
# by identifying the mode of all month values per TC best track
most_frequent_month = all_TC_tracks.groupby('SID')['MONTH'].agg(lambda x: mode(x)[0])

most_frequent_month = most_frequent_month.map({
    1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May',
    6: 'June', 7: 'July', 8: 'August', 9: 'September', 10: 'October',
    11: 'November', 12: 'December'})

all_TC_tracks = all_TC_tracks.merge(
    most_frequent_month, how='left', on='SID').rename(
        columns={'MONTH_y': 'MOST_FREQUENT_MONTH'})

all_TC_tracks.drop(columns=['ISO_TIME','MONTH_x'], inplace=True)
all_TC_tracks.rename(columns={'MOST_FREQUENT_MONTH': 'MONTH'}, inplace=True)

def plot_monthly_tracks(month):
    # Filter the DataFrame to only include rows for a specific month
    monthly_tracks = all_TC_tracks[all_TC_tracks['MONTH'] == month]

    # Group the DataFrame by 'SID'
    grouped = monthly_tracks.groupby('SID')

    # Initialize an empty list to store the dictionaries
    monthly_TCs = []
    
    # Iterate over the grouped DataFrame
    for name, group in grouped:
        # Create a dictionary where the key is the SID and the value is a list of lat-lon tuples
        tc_dict = {name: list(zip(group['LAT'], group['LON']))}
        # Append the dictionary to the list
        monthly_TCs.append(tc_dict)
    
    # Create a 2D array for the histogram
    hist = np.zeros((50, 75)) #50S-0N, 20-95E

    # Iterate over the tropical cyclones
    for tc in monthly_TCs:
        for sid, locations in tc.items():
            # Create a set to store the grid boxes that this tropical cyclone passes through
            visited_boxes = set()
            for lat, lon in locations:
                # Calculate the indices of the grid box that this location falls into
                lat_index = int(lat - (-50))
                lon_index = int(lon - 20)
                # Check if the indices are within the expected range
                if 0 <= lat_index < 50 and 0 <= lon_index < 75:
                    # Add the grid box to the set of visited boxes
                    visited_boxes.add((lat_index, lon_index))
            # Increment the count in the histogram for each visited box
            for lat_index, lon_index in visited_boxes:
                hist[lat_index, lon_index] += 1

    # Define the edges of the bins for the histogram
    xedges = np.arange(20, 96, 1)  # 75 elements
    yedges = np.arange(-50, 1, 1)  # 50 elements

    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([20, 95, -50, 0], ccrs.PlateCarree())

    gridlines = ax.gridlines(draw_labels=False, xlocs=np.arange(-180, 181, 5), ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
    xticks = np.arange(20, 96, 5)
    yticks = np.arange(-50, 1, 5)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    ax.coastlines('10m', edgecolor='black', linewidth=2.5)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2.5)

    # Plot the histogram as a heatmap on the geographical map
    cax = ax.pcolormesh(xedges, yedges, hist, vmin=0, vmax=50, cmap='PuRd', shading='auto')

    # Add a colorbar
    cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    norm = Normalize(vmin=0, vmax=50)
    colorbar = ColorbarBase(cax, cmap='PuRd', norm=norm, extend='max', ticks=np.arange(0, 51, 2))
    colorbar.set_label('Frequency of tropical cyclone passage', size=25)
    cax.tick_params(labelsize=15)

    # Iterate over the tropical cyclones
    for tc in monthly_TCs:
        for sid, locations in tc.items():
            # Unzip the list of lat-lon tuples into two separate lists
            lats, lons = zip(*locations)
            # Plot the track of the tropical cyclone
            ax.plot(lons, lats, color='black', linewidth=0.2, transform=ccrs.PlateCarree())

    # Add a title and subtitle
    fig.suptitle('IBTrACS | Monthly climatology of Southwest Indian Ocean tropical cyclone tracks (1848–2023)', fontsize=28, y=0.803)
    ax.set_title(month.upper(), fontsize=40, pad=20)

    plt.show()

# List of months
months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

# Generate a plot for each month
for month in months:
    plot_monthly_tracks(month)
    # break

# End of program