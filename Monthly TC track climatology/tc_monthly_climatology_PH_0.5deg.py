# Program for generating monthly climatology maps of Philippine tropical cyclone tracks
# using IBTrACS data

import numpy as np
import pandas as pd
from scipy.stats import mode
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Extracting TC best track data from the IBTrACS CSV file
all_TC_tracks = pd.read_csv('C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/ibtracs.WP.list.v04r00.csv',
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

# Define the Philippine Area of Responsibility
PAR_corners_coords = [(120, 25), (135, 25), (135, 5), (115, 5), (115, 15), (120, 21)]
PAR = Polygon(PAR_corners_coords)

# Identifying if a TC center lat-lon location is within the PAR
def is_in_PAR(lat, lon):
    point = Point(lon, lat)  # Note that Point takes (lon, lat) not (lat, lon)
    return PAR.contains(point)

all_TC_tracks['IN_REGION'] = all_TC_tracks.apply(lambda row: is_in_PAR(row['LAT'], row['LON']), axis=1)

# Removing TCs that did not track within the PAR
sids_in_region = all_TC_tracks.loc[all_TC_tracks['IN_REGION'], 'SID'].unique()
all_TC_tracks = all_TC_tracks[all_TC_tracks['SID'].isin(sids_in_region)]
all_TC_tracks.drop(columns=['IN_REGION'], inplace=True)

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
    hist = np.zeros((60, 88))

    # Iterate over the tropical cyclones
    for tc in monthly_TCs:
        for sid, locations in tc.items():
            # Create a set to store the grid boxes that this tropical cyclone passes through
            visited_boxes = set()
            for lat, lon in locations:
                # Calculate the indices of the grid box that this location falls into
                # Multiply by 2 to account for the 0.5 degree grid size
                lat_index = int((lat - 0) * 2)
                lon_index = int((lon - 106) * 2)
                # Check if the indices are within the expected range
                if 0 <= lat_index < 60 and 0 <= lon_index < 88:
                    # Add the grid box to the set of visited boxes
                    visited_boxes.add((lat_index, lon_index))
            # Increment the count in the histogram for each visited box
            for lat_index, lon_index in visited_boxes:
                hist[lat_index, lon_index] += 1

    # Define the edges of the bins for the histogram
    xedges = np.arange(106, 150.5, 0.5)  # 88 elements
    yedges = np.arange(0, 30.5, 0.5)  # 60 elements

    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([106, 150, 0, 30], ccrs.PlateCarree())

    gridlines = ax.gridlines(draw_labels=False, xlocs=np.arange(-180, 181, 1), ylocs=np.arange(-90, 91, 1), color='gray', linestyle='--')
    xticks = np.arange(106, 151, 2)
    yticks = np.arange(0, 31, 2)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    ax.coastlines('10m', edgecolor='black', linewidth=2.5)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2.5)

    lons = [120, 135, 135, 115, 115, 120, 120]
    lats = [25, 25, 5, 5, 15, 21, 25]
    ax.plot(lons, lats, color='black', linewidth=5, transform=ccrs.PlateCarree(), linestyle='--')

    # Plot the histogram as a heatmap on the geographical map
    cax = ax.pcolormesh(xedges, yedges, hist, vmin=0, vmax=25, cmap='PuRd', shading='auto')

    # Add a colorbar
    cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    norm = Normalize(vmin=0, vmax=25)
    colorbar = ColorbarBase(cax, cmap='PuRd', norm=norm, extend='max', ticks=np.arange(0, 26, 1))
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
    fig.suptitle('IBTrACS | Monthly climatology of Philippine tropical cyclone tracks (1884–2023)', fontsize=30, y=0.815)
    ax.set_title(month.upper(), fontsize=40, pad=20)

    plt.show()

# List of months
months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

# Generate a plot for each month
for month in months:
    plot_monthly_tracks(month)
    # break

# End of program